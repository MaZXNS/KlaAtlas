"""KlaAtlas analysis routines."""
from __future__ import annotations

import argparse,json

from pathlib import Path

import numpy as np,pandas as pd

from scipy.stats import norm,t

from statsmodels.stats.meta_analysis import combine_effects

from statsmodels.stats.multitest import multipletests

AA='ACDEFGHIKLMNPQRSTVWY'

def prepare(f):
    keys=['protein_unit_id','position_decile'];f=f.sort_values(keys).copy();grp=f.groupby(keys,sort=False);code=grp.ngroup().to_numpy();meta=grp[['taxon_id','protein_unit_id']].first();n=np.bincount(code).astype(float);m=np.bincount(code,weights=f.is_kla).astype(float);protein=pd.factorize(meta.protein_unit_id,sort=False)[0];taxon=meta.taxon_id.astype(str).to_numpy();return f.reset_index(drop=True),code,n,m,protein,taxon

def run(f,name):
    f,g,n,m,protein,taxon=prepare(f);y=f.is_kla.to_numpy(bool);E=np.frombuffer(''.join(f.window).encode(),dtype=np.uint8).reshape(len(f),21);eligible=[x for x in sorted(set(taxon)) if len(np.unique(protein[taxon==x]))>=200 and m[taxon==x].sum()>=200];taxrows=[];summ=[]
    for dist in range(1,11):
        for aa in AA:
            sides=[]
            for col in [10-dist,10+dist]:
                tag=E[:,col]==ord(aa);tot=np.bincount(g,weights=tag,minlength=len(n));A=np.bincount(g,weights=tag&y,minlength=len(n));B=m-A;C=tot-A;D=n-m-tot+A;sides.append((A*D/n,B*C/n))
            local=[]
            for tax in eligible:
                ix=taxon==tax;ests=[];scores=[]
                for num,den in sides:
                    nu=num[ix];de=den[ix]
                    if nu.sum()<=0 or de.sum()<=0:break
                    ests.append(np.log2(nu.sum()/de.sum()));scores.append(np.bincount(protein[ix],weights=(nu/nu.sum()-de/de.sum())/np.log(2),minlength=protein.max()+1))
                row={'population':name,'taxon_id':tax,'distance':dist,'residue':aa,'estimable':False}
                if len(ests)==2:
                    N=len(np.unique(protein[ix]));score=scores[1]-scores[0];se=np.sqrt((score**2).sum()*N/(N-1));delta=ests[1]-ests[0];row.update(estimable=bool(se>0),right_minus_left=delta,standard_error=se,ci_low=delta-1.96*se,ci_high=delta+1.96*se,p_value=2*norm.sf(abs(delta/se)) if se>0 else np.nan,proteins=N,Kla=int(m[ix].sum()),K=int(n[ix].sum()));
                    if se>0:local.append(row)
                taxrows.append(row)
            summary={'population':name,'distance':dist,'residue':aa,'taxa':len(local),'estimable':False,'reason':'fewer_than_3_taxa'}
            if len(local)>=3:
                eff=np.array([r['right_minus_left'] for r in local]);var=np.array([r['standard_error']**2 for r in local]);fit=combine_effects(eff,var,method_re='iterated',use_t=True);delta=float(fit.mean_effect_re);se=float(fit.sd_eff_w_re_hksj);df=len(local)-1;crit=t.ppf(.975,df);summary.update(estimable=bool(np.isfinite(se) and se>0),reason='estimable',right_minus_left=delta,standard_error=se,ci_low=delta-crit*se,ci_high=delta+crit*se,p_value=2*t.sf(abs(delta/se),df) if se>0 else np.nan,i_squared=max(0,float(fit.i2)))
            summ.append(summary)
        print(name,dist,flush=True)
    S=pd.DataFrame(summ)
    if 'p_value' not in S:S['p_value']=np.nan
    S['q_BH200']=multipletests(S.p_value.fillna(1),method='fdr_bh')[1];S.loc[~S.estimable,'q_BH200']=np.nan;return pd.DataFrame(taxrows),S,{'population':name,'K':len(f),'Kla':int(f.is_kla.sum()),'proteins':f.protein_unit_id.nunique(),'eligible_taxa':eligible,'family':200}
