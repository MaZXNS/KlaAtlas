"""KlaAtlas analysis routines."""
from __future__ import annotations

import argparse

import ctypes

import hashlib

from pathlib import Path

import subprocess

import sys

import numpy as np

import pandas as pd

class MultiscaleKernel:
    """Compile the fixed-count permutation kernel."""
    def __init__(self,cache:Path):
        source=Path(__file__).with_name('multiscale_kernel.c')
        self.sha=hashlib.sha256(source.read_bytes()).hexdigest()
        cache.mkdir(parents=True,exist_ok=True)
        path=cache/f'multiscale_{self.sha[:16]}{".dylib" if sys.platform=="darwin" else ".so"}'
        if not path.exists():subprocess.run(['cc','-O3','-std=c99','-shared','-fPIC',str(source),'-lm','-o',str(path)],check=True)
        self.lib=ctypes.CDLL(str(path.resolve()))
        i32=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        f64=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        self.lib.multiscale_eval.argtypes=[i32,ctypes.c_int,i32,ctypes.c_int,i32,ctypes.c_int,ctypes.c_int,ctypes.c_uint64,f64]
        self.lib.multiscale_eval.restype=ctypes.c_int

    def evaluate(self,k:np.ndarray,kla:np.ndarray,widths:list[int],replicates:int,seed:int) -> tuple[dict,np.ndarray]:
        """Return one protein's statistics and centered exchangeable null draws."""
        k=np.ascontiguousarray(np.sort(k),dtype=np.int32);kla=np.ascontiguousarray(np.sort(kla),dtype=np.int32)
        w=np.ascontiguousarray(widths,dtype=np.int32)
        if len(np.unique(k))!=len(k) or not set(kla)<=set(k):raise ValueError('Invalid K coordinate opportunity set')
        base=6+3*len(w);result=np.zeros(base+replicates*(len(w)+1))
        status=self.lib.multiscale_eval(k,len(k),kla,len(kla),w,len(w),replicates,seed,result)
        if status:raise RuntimeError(f'Multiscale permutation failure {status}')
        effect=float(np.log2(result[0])-result[1])
        record={'observed_mean_nearest':result[0],'null_mean_log2_nearest':result[1],
            'centered_log2_nearest':effect,'geometric_normalized_nearest_ratio':2**effect,
            'p_nearest':result[3],'arithmetic_null_mean_nearest':result[4],'median_null_nearest':result[5]}
        for j,width in enumerate(widths):
            record.update({f'observed_scan_{width}':result[6+j],f'expected_scan_{width}':result[6+len(w)+j],
                f'fraction_excess_{width}':(result[6+j]-result[6+len(w)+j])/len(kla),f'p_scan_{width}':result[6+2*len(w)+j]})
        record['p_scan_family']=min(1.,len(w)*min(record[f'p_scan_{width}'] for width in widths))
        return record,result[base:].reshape(replicates,len(w)+1)
