"""KlaAtlas analysis routines."""
import exact_conditional as exact
from pathlib import Path

import sys

import numpy as np,pandas as pd

import exact_conditional as exact

BASE=['kr_fraction','de_fraction','hydrophobicity','within_decile_position','theory_available','theory_length_log1p']

def cluster_meat(scores:np.ndarray,labels:np.ndarray)->tuple[np.ndarray,int]:
 """Aggregate stratum scores by an explicitly defined whole cluster."""
 codes,unique=pd.factorize(labels);summed=np.zeros((len(unique),scores.shape[1]));np.add.at(summed,codes,scores);V=summed.T@summed
 if len(unique)>1:V*=len(unique)/(len(unique)-1)
 return V,len(unique)

def fit_two_way(frame:pd.DataFrame,features:list[str],kernel:object,protein_clusters:dict[str,str],maxiter:int=300)->tuple[dict,object,dict]:
 """Exact source×protein×position conditioning; protein and reuse-component sandwich."""
 data=exact.prepare_data(frame,features);fit=exact.fit_conditional(kernel,data,maxiter=maxiter);fake=data.protein_ids[data.stratum_protein];comp=np.array([x.split('@@',1)[0] for x in fake]);original=np.array([x.split('@@',1)[1] for x in fake]);protein=np.array([protein_clusters.get(x,x) for x in original]);vp,np_=cluster_meat(fit['stratum_scores'],protein);vs,ns=cluster_meat(fit['stratum_scores'],comp);vi,ni=cluster_meat(fit['stratum_scores'],np.char.add(np.char.add(comp,'@@'),protein));H=fit['hessian_inverse'];cov=H@(vp+vs-vi)@H;ok=bool(fit['success'] and fit['gradient_per_stratum_max']<1e-5 and fit['hessian_min_eigenvalue']>1e-7 and abs(fit['beta']).max()<20)
 meta={'status':'estimable' if ok else 'failed_diagnostic','features':data.features,'beta':fit['beta'].tolist(),'covariance_two_way':cov.tolist(),'covariance_min_eigenvalue':float(np.linalg.eigvalsh((cov+cov.T)/2).min()),'K_rows':len(data.y),'case_rows':int(data.y.sum()),'source_components':ns,'sequence_protein_clusters':np_,'intersection_clusters':ni,'strata':len(data.starts)-1,'loss':fit['loss'],'gradient_per_stratum':fit['gradient_per_stratum_max'],'hessian_min':fit['hessian_min_eigenvalue']};fit['covariance_two_way']=cov;return meta,data,fit
