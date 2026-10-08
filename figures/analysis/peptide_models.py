"""KlaAtlas analysis routines."""
from __future__ import annotations

import argparse,json,hashlib,shutil,sys,importlib.util,itertools

from pathlib import Path

import numpy as np

import pandas as pd

import patsy

from scipy.stats import norm,chi2

from statsmodels.stats.multitest import multipletests

BASE=['kr_fraction','de_fraction','hydrophobicity','within_decile_position']

KEYS=['protein_unit_id','reference_sequence_id','k_position']

def save(d:pd.DataFrame,p:Path)->None:
    d.to_csv(p,sep='\t',index=False,compression={'method':'gzip','mtime':0} if p.suffix=='.gz' else None)

def info(f:pd.DataFrame)->pd.DataFrame:
    counts=f.groupby(['protein_unit_id','position_decile']).is_kla.agg(['size','sum']);valid=counts[(counts['sum']>0)&(counts['sum']<counts['size'])].index;return f[pd.MultiIndex.from_frame(f[['protein_unit_id','position_decile']]).isin(valid)].copy().reset_index(drop=True)

def fitting(frame:pd.DataFrame,features:list,name:str,out:Path,exact:object,kernel:object)->dict:
    """Fit exact multi-case conditional likelihood; retain complete covariance and cluster influence."""
    f=info(frame);meta=dict(model=name,opportunity_rows=len(frame),informative_rows=len(f),reported_rows=int(f.is_kla.sum()),unique_K=f.coordinate_id.nunique(),unique_cases=f.loc[f.is_kla.eq(1),'coordinate_id'].nunique(),proteins=f.protein_unit_id.nunique(),strata=f.groupby(['protein_unit_id','position_decile']).ngroups,requested_features=features);rows=[];result={'meta':meta}
    try:
        d=exact.prepare_data(f,features);z=exact.fit_conditional(kernel,d,maxiter=350);rank=int(np.linalg.matrix_rank(d.x-pd.DataFrame(d.x).groupby(np.repeat(np.arange(len(d.starts)-1),np.diff(d.starts))).transform('mean').to_numpy()));ok=z['success'] and z['gradient_per_stratum_max']<1e-5 and z['hessian_min_eigenvalue']>1e-7 and np.max(abs(z['beta']))<20 and rank==len(d.features)
        meta.update(status='estimable' if ok else 'failed_rank_or_fit_diagnostic',active_features=d.features,beta=z['beta'].tolist(),covariance_cluster=z['covariance_cluster'].tolist(),gradient_per_stratum=z['gradient_per_stratum_max'],hessian_min_eigenvalue=z['hessian_min_eigenvalue'],within_stratum_rank=rank,loglik=-z['loss'])
        result.update(beta=z['beta'],cov=z['covariance_cluster'],features=d.features,proteins=d.protein_ids,influence=z['protein_scores']@z['hessian_inverse'].T,valid=bool(ok))
        for term in features:
            r=dict(model=name,feature=term,estimable=bool(ok and term in d.features),reason=meta['status'] if term in d.features else 'invariant')
            if r['estimable']:
                j=d.features.index(term);b=z['beta'][j];se=np.sqrt(z['covariance_cluster'][j,j]);r.update(log2_or=b/np.log(2),standard_error_log2=se/np.log(2),ci95_lower=(b-1.96*se)/np.log(2),ci95_upper=(b+1.96*se)/np.log(2),p_value=2*norm.sf(abs(b/se)),p_value_log10=(np.log(2)+norm.logsf(abs(b/se)))/np.log(10))
            rows.append(r)
        if ok:save(pd.DataFrame(result['influence'],columns=d.features).assign(protein_unit_id=d.protein_ids),out/(name+'_PROTEIN_INFLUENCE.tsv.gz'))
    except Exception as ex:
        meta.update(status='failed',reason=str(ex));result['valid']=False;rows=[dict(model=name,feature=t,estimable=False,reason=str(ex)) for t in features]
    save(pd.DataFrame(rows),out/(name+'_EFFECTS.tsv'));save(f,out/(name+'_INFORMATIVE_ROWS.tsv.gz'));(out/(name+'_MODEL.json')).write_text(json.dumps(meta,indent=2,default=str));print(name,meta['status'],meta['reported_rows'],meta['proteins'],flush=True);return result

def spline_frames(frames:dict,column:str,prefix:str,knots:list)->tuple[dict,object]:
    """Create a centered natural-cubic basis from the opportunity covariates."""
    x=np.concatenate([np.log1p(f[column].to_numpy()) for f in frames.values()]);lo=float(x.min());hi=float(x.max());kn=tuple(float(np.log1p(z)) for z in knots if lo<np.log1p(z)<hi);dm=patsy.dmatrix('cr(x, knots=kn, lower_bound=lo, upper_bound=hi, constraints="center") - 1',dict(x=x,kn=kn,lo=lo,hi=hi));design=dm.design_info;start=0;names=[prefix+str(i) for i in range(dm.shape[1])];out={}
    for label,f in frames.items():
        g=f.copy();arr=np.asarray(dm)[start:start+len(f)];start+=len(f)
        for i,n in enumerate(names):g[n]=arr[:,i]
        out[label]=g
    return out,(design,names,dict(lo=lo,hi=hi,kn=kn))

def curve(model:dict,basis:tuple,grid:np.ndarray,reference:float,name:str)->pd.DataFrame:
    """Cluster-Wald relative odds curve within observed covariate support."""
    design,names,params=basis;xx=np.log1p(np.r_[grid,reference]);mat=np.asarray(patsy.build_design_matrices([design],dict(x=xx,**params))[0]);D=mat[:-1]-mat[-1];rows=[]
    for x,v in zip(grid,D):
        row=dict(model=name,x=x,reference=reference,estimable=model.get('valid',False))
        if row['estimable'] and all(n in model['features'] for n in names):
            full=np.zeros(len(model['features']))
            for val,n in zip(v,names):full[model['features'].index(n)]=val
            be=float(full@model['beta']);se=float(np.sqrt(max(0,full@model['cov']@full)));row.update(log2_relative_or=be/np.log(2),ci95_lower=(be-1.96*se)/np.log(2),ci95_upper=(be+1.96*se)/np.log(2))
        else:row['estimable']=False
        rows.append(row)
    return pd.DataFrame(rows)
