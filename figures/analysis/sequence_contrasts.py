"""KlaAtlas analysis routines."""
from __future__ import annotations

import argparse

from pathlib import Path

import numpy as np

import pandas as pd

from scipy.stats import norm,t as student_t
from statsmodels.stats.multitest import multipletests
def bh_adjust(p, family):
    p=np.asarray(p,float); full=np.ones(family); full[:len(p)]=np.nan_to_num(p,nan=1); return multipletests(full,method='fdr_bh')[1][:len(p)]

from statsmodels.stats.meta_analysis import combine_effects

AA='ACDEFGHIKLMNPQRSTVWY'

def prepare(frame: pd.DataFrame) -> tuple[pd.DataFrame,np.ndarray,np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    """Enumerate informative matched strata and protein clusters without control reuse."""
    keys=['protein_unit_id','position_decile']
    frame=frame.sort_values(keys).copy()
    groups=frame.groupby(keys,sort=False).is_kla.agg(['size','sum'])
    valid=groups[(groups['sum']>0)&(groups['sum']<groups['size'])].reset_index()[keys]
    frame=frame.merge(valid,on=keys,how='inner',validate='many_to_one').sort_values(keys).reset_index(drop=True)
    grouped=frame.groupby(keys,sort=False)
    code=grouped.ngroup().to_numpy()
    groupmeta=grouped[['taxon_id','protein_unit_id']].first()
    protein=pd.factorize(groupmeta.protein_unit_id,sort=False)[0]
    taxon=groupmeta.taxon_id.to_numpy(str)
    n=np.bincount(code).astype(float);m=np.bincount(code,weights=frame.is_kla).astype(float)
    return frame,code,n,m,protein,taxon

def sequence_effects(frame: pd.DataFrame, scenario: str) -> tuple[pd.DataFrame,pd.DataFrame,dict]:
    """Estimate protein-decile MH effects with protein-cluster uncertainty and taxon meta-analysis."""
    frame,g,n,m,protein,taxon=prepare(frame)
    y=frame.is_kla.to_numpy(dtype=bool)
    encoded=np.frombuffer(''.join(frame.window).encode('ascii'),dtype=np.uint8).reshape(len(frame),21)
    eligible_taxa=[]
    for tax in sorted(set(taxon)):
        sel=taxon==tax
        if len(np.unique(protein[sel]))>=200 and m[sel].sum()>=200: eligible_taxa.append(tax)
    tax_rows=[];meta_rows=[]
    for offset in list(range(-10,0))+list(range(1,11)):
        col=encoded[:,offset+10]
        for residue in AA:
            annotated=col==ord(residue)
            total=np.bincount(g,weights=annotated,minlength=len(n))
            a=np.bincount(g,weights=annotated&y,minlength=len(n))
            b=m-a;c=total-a;d=n-m-total+a
            numerator=a*d/n;denominator=b*c/n
            local=[]
            for tax in eligible_taxa:
                select=taxon==tax;num=numerator[select];den=denominator[select]
                if num.sum()<=0 or den.sum()<=0: continue
                beta=np.log2(num.sum()/den.sum())
                influence=(num/num.sum()-den/den.sum())/np.log(2)
                score=np.bincount(protein[select],weights=influence)
                proteins=len(np.unique(protein[select]));se=np.sqrt((score**2).sum()*proteins/(proteins-1))
                row={'scenario':scenario,'taxon_id':tax,'offset':offset,'residue':residue,
                    'log2_or':beta,'standard_error':se,'ci95_lower':beta-1.96*se,'ci95_upper':beta+1.96*se,
                    'p_value':2*norm.sf(abs(beta/se)) if se>0 else np.nan,
                    'proteins':proteins,'kla_sites':int(m[select].sum()),'other_k':int((n[select]-m[select]).sum())}
                tax_rows.append(row);local.append(row)
            if len(local)>=3:
                effects=np.array([r['log2_or'] for r in local]);variances=np.array([r['standard_error']**2 for r in local])
                combined=combine_effects(effects,variances,method_re='iterated',use_t=True)
                beta=float(combined.mean_effect_re);se=float(combined.sd_eff_w_re_hksj);df=len(local)-1
                if se>0 and np.isfinite(se):
                    crit=student_t.ppf(.975,df)
                    meta_rows.append({'scenario':scenario,'offset':offset,'residue':residue,'taxa':len(local),
                        'log2_or':beta,'standard_error':se,'ci95_lower':beta-crit*se,'ci95_upper':beta+crit*se,
                        'p_value':2*student_t.sf(abs(beta/se),df),'i_squared':max(0,float(combined.i2)),
                        'direction_concordance':max(float((effects>0).mean()),float((effects<0).mean()))})
        print(f'{scenario}: offset {offset}',flush=True)
    taxout=pd.DataFrame(tax_rows);metaout=pd.DataFrame(meta_rows)
    taxout['q_bh_400']=np.nan
    for tax,indices in taxout.groupby('taxon_id').groups.items():
        taxout.loc[indices,'q_bh_400']=bh_adjust(taxout.loc[indices,'p_value'].to_numpy(),400)
    metaout['q_bh_400']=bh_adjust(metaout.p_value.to_numpy(),400)
    return taxout,metaout,{'scenario':scenario,'k_count':len(frame),'kla_count':int(y.sum()),
        'other_k':int((~y).sum()),'proteins':frame.protein_unit_id.nunique(),'strata':len(n),
        'estimable_taxa':eligible_taxa,'family_size':400}
