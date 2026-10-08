"""KlaAtlas analysis routines."""
from __future__ import annotations

import ctypes

import hashlib

import subprocess

import sys

from dataclasses import dataclass

from pathlib import Path

from typing import Any

import numpy as np

import pandas as pd

from scipy.optimize import minimize

@dataclass
class ConditionalData:
    """Contiguous K rows and boundaries of informative matched strata."""
    x: np.ndarray
    y: np.ndarray
    starts: np.ndarray
    stratum_protein: np.ndarray
    protein_ids: np.ndarray
    features: list[str]
    row_index: np.ndarray

class ExactConditional:
    """Stable exact conditional likelihood, accelerated with a local C kernel."""

    def __init__(self, cache_dir: Path):
        """Compile the conditional-likelihood kernel."""
        source=Path(__file__).with_name('conditional_kernel.c')
        digest=hashlib.sha256(source.read_bytes()).hexdigest()[:16]
        cache_dir.mkdir(parents=True,exist_ok=True)
        extension='.dylib' if sys.platform=='darwin' else '.so'
        library=cache_dir/f'conditional_{digest}{extension}'
        if not library.exists():
            subprocess.run(['cc','-O3','-std=c99','-shared','-fPIC',str(source),'-lm','-o',str(library)],check=True)
        self.library_path=library
        self.lib=ctypes.CDLL(str(library.resolve()))
        f64=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        i32=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        self.lib.conditional_eval.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_int,
            f64,i32,i32,f64,f64,f64,f64,f64,f64]
        self.lib.conditional_eval.restype=ctypes.c_int

    def evaluate(self, data: ConditionalData, beta: np.ndarray,
                 weights: np.ndarray | None=None) -> tuple[float,np.ndarray,np.ndarray,np.ndarray]:
        """Return total negative log likelihood, gradient, stratum scores and log likelihoods."""
        n,p=data.x.shape; g=len(data.starts)-1
        beta=np.ascontiguousarray(beta,dtype=np.float64)
        if weights is None: weights=np.ones(g,dtype=np.float64)
        weights=np.ascontiguousarray(weights,dtype=np.float64)
        loss=np.zeros(1); gradient=np.zeros(p); scores=np.zeros((g,p)); ll=np.zeros(g)
        status=self.lib.conditional_eval(n,p,g,data.x,data.y,data.starts,beta,weights,
                                        loss,gradient,scores,ll)
        if status:
            raise FloatingPointError(f'Conditional recursion failed with status {status}')
        return float(loss[0]),gradient,scores,ll

def prepare_data(frame: pd.DataFrame, features: list[str]) -> ConditionalData:
    """Prepare matched strata and variable features."""
    ordered=frame.sort_values(['protein_unit_id','position_decile','k_position']).copy()
    ordered['_input_row_index']=ordered.index
    keys=['protein_unit_id','position_decile']
    counts=ordered.groupby(keys,sort=False).is_kla.agg(['size','sum'])
    valid=counts[(counts['sum']>0)&(counts['sum']<counts['size'])].reset_index()[keys]
    ordered=ordered.merge(valid,on=keys,how='inner',validate='many_to_one')
    groups=ordered.groupby(keys,sort=False)
    lengths=groups.size().to_numpy(dtype=np.int32)
    starts=np.r_[0,np.cumsum(lengths)].astype(np.int32)
    means=groups[features].transform('mean')
    varied=((ordered[features]-means).pow(2).sum()>1e-12)
    active=[feature for feature in features if varied[feature]]
    if not active or not len(lengths):
        raise ValueError('No identifiable feature/informative strata')
    stratum_names=np.asarray([key[0] for key in groups.groups],dtype=str)
    protein_ids,protein_codes=np.unique(stratum_names,return_inverse=True)
    return ConditionalData(
        np.ascontiguousarray(ordered[active].to_numpy(dtype=float)),
        np.ascontiguousarray(ordered.is_kla.to_numpy(dtype=np.int32)),starts,
        protein_codes,protein_ids,active,ordered['_input_row_index'].to_numpy())

def fit_conditional(kernel: ExactConditional, data: ConditionalData,
                    start: np.ndarray | None=None, weights: np.ndarray | None=None,
                    maxiter: int=250) -> dict[str,Any]:
    """Fit the exact likelihood and obtain a numerical Hessian and cluster scores."""
    p=data.x.shape[1]; g=len(data.starts)-1
    norm=float(g if weights is None else np.sum(weights))
    if start is None: start=np.zeros(p)
    scale=np.maximum(data.x.std(axis=0),1e-6)
    def objective(theta: np.ndarray) -> tuple[float,np.ndarray]:
        beta=theta/scale
        loss,gradient,_,_=kernel.evaluate(data,beta,weights)
        return loss/norm,gradient/(norm*scale)
    result=minimize(objective,start*scale,jac=True,method='L-BFGS-B',
                    options={'maxiter':maxiter,'ftol':1e-12,'gtol':1e-8,'maxls':40})
    beta=result.x/scale
    loss,gradient,scores,ll=kernel.evaluate(data,beta,weights)
    hessian=np.empty((p,p))
    for j in range(p):
        step=1e-4*max(1,abs(beta[j])); delta=np.zeros(p); delta[j]=step
        gp=kernel.evaluate(data,beta+delta,weights)[1]
        gm=kernel.evaluate(data,beta-delta,weights)[1]
        hessian[:,j]=(gp-gm)/(2*step)
    hessian=(hessian+hessian.T)/2
    eigenvalues=np.linalg.eigvalsh(hessian)
    inverse=np.linalg.pinv(hessian,rcond=1e-10)
    protein_scores=np.zeros((len(data.protein_ids),p))
    np.add.at(protein_scores,data.stratum_protein,scores)
    covariance=inverse@(protein_scores.T@protein_scores)@inverse
    if len(protein_scores)>1: covariance*=len(protein_scores)/(len(protein_scores)-1)
    return {'beta':beta,'loss':loss,'gradient':gradient,'hessian':hessian,'hessian_inverse':inverse,
            'covariance_cluster':covariance,'protein_scores':protein_scores,'stratum_scores':scores,
            'stratum_loglik':ll,'success':bool(result.success),
            'message':str(result.message),'iterations':result.nit,
            'gradient_per_stratum_max':float(np.abs(gradient).max()/g),
            'hessian_min_eigenvalue':float(eigenvalues.min()),'condition_number':float(np.linalg.cond(hessian))}

def protein_bootstrap(kernel: ExactConditional, data: ConditionalData,
                      fit: dict[str,Any], replicates: int, seed: int,
                      full_refit: bool=True) -> tuple[np.ndarray,list[dict[str,Any]]]:
    """Resample whole proteins, keeping every decile from a sampled protein together."""
    rng=np.random.default_rng(seed); n=len(data.protein_ids)
    estimates=[]; diagnostics=[]
    for b in range(replicates):
        counts=np.bincount(rng.integers(0,n,n),minlength=n).astype(float)
        if full_refit:
            weights=counts[data.stratum_protein]
            normalization=weights.sum()
            def objective(beta: np.ndarray) -> tuple[float,np.ndarray]:
                loss,gradient,_,_=kernel.evaluate(data,beta,weights)
                return loss/normalization,gradient/normalization
            result=minimize(objective,fit['beta'],jac=True,method='L-BFGS-B',
                options={'maxiter':100,'ftol':1e-10,'gtol':2e-7,'maxls':30})
            estimate=result.x
            diagnostics.append({'replicate':b,'success':bool(result.success)})
        else:
            estimate=fit['beta']+fit['hessian_inverse']@((counts-1)@fit['protein_scores'])
            diagnostics.append({'replicate':b,'success':True,'method':'one_step_cluster_score_approximation'})
        estimates.append(estimate)
    return np.asarray(estimates),diagnostics
