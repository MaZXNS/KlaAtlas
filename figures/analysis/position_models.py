"""KlaAtlas analysis routines."""
from __future__ import annotations

import argparse, hashlib, json

from pathlib import Path

import numpy as np

import pandas as pd

from statsmodels.stats.multitest import multipletests

SCOPES={
    "global_unknown":("global_is_kla","global_is_censored"),
}

def sha(p:Path)->str:
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def save(d:pd.DataFrame,p:Path)->None:
    d.to_csv(p,sep='\t',index=False,compression={'method':'gzip','mtime':0} if p.suffix=='.gz' else None)

def analyze(k:pd.DataFrame,scope:str,case_col:str,censor_col:str,out:Path,draws:int,seed:int)->dict:
    """Build protein-balanced 20-bin profiles and fixed-burden nulls."""
    rows=[];den=[];excluded=[]
    for pid,g0 in k.groupby('protein_unit_id',sort=True):
        g=g0.loc[~g0[censor_col].eq(1)].copy();case=g[case_col].eq(1).to_numpy();n=len(g);m=int(case.sum())
        if m==0 or m>n:
            excluded.append({'protein_unit_id':pid,'reason':'no_scope_case' if m==0 else 'case_denominator_error','K':n,'Kla':m});continue
        rel=(g.k_position.to_numpy(int)-1)/(g.sequence_length.iloc[0]-1);bi=np.minimum((rel*20).astype(int),19)
        kc=np.bincount(bi,minlength=20);cc=np.bincount(bi[case],minlength=20);expect=kc/n;obs=cc/m
        den.append({'protein_unit_id':pid,'taxon_id':str(g.taxon_id.iloc[0]),'K':n,'Kla':m,'null_variable':bool(m<n and np.count_nonzero(kc)>1)})
        rows.extend({'protein_unit_id':pid,'taxon_id':str(g.taxon_id.iloc[0]),'bin':j,'K_bin':int(kc[j]),'Kla_bin':int(cc[j]),'observed_fraction':obs[j],'expected_fraction':expect[j],'difference':obs[j]-expect[j]} for j in range(20))
    d=pd.DataFrame(den);b=pd.DataFrame(rows);save(d,out/f'{scope}_protein_denominators.tsv');save(b,out/f'{scope}_protein_bins.tsv.gz');save(pd.DataFrame(excluded),out/f'{scope}_exclusions.tsv')
    groups={'all':np.ones(len(d),bool),**{f'taxon:{t}':d.taxon_id.eq(t).to_numpy() for t in sorted(d.taxon_id.unique())}}
    null=np.zeros((len(groups),draws,20));names=list(groups);lookup={x:i for i,x in enumerate(names)};rng=np.random.default_rng(seed)
    nmat=b.K_bin.to_numpy().reshape(-1,20);cmat=b.Kla_bin.to_numpy().reshape(-1,20)
    for i,r in enumerate(d.itertuples(index=False)):
        z=rng.multivariate_hypergeometric(nmat[i],r.Kla,size=draws)/r.Kla-nmat[i]/r.K
        null[lookup['all']]+=z;null[lookup[f'taxon:{r.taxon_id}']]+=z
        if i and i%5000==0:print(scope,'P042 proteins',i,flush=True)
    summaries=[];tests=[];delta=b.difference.to_numpy().reshape(-1,20)
    for name,mask in groups.items():
        n=int(mask.sum());z=null[lookup[name]]/max(n,1);diff=delta[mask].mean(0);exp=(nmat[mask]/d.K.to_numpy()[mask,None]).mean(0);lo,hi=np.quantile(z,[.025,.975],axis=0)
        stat=float(np.square(diff).sum());p=(1+(np.square(z).sum(1)>=stat).sum())/(draws+1)
        tests.append({'scope':scope,'group':name,'proteins':n,'K':int(d.loc[mask,'K'].sum()),'Kla':int(d.loc[mask,'Kla'].sum()),'omnibus_statistic':stat,'p_value':p,'terminal_difference':float(diff[[0,19]].sum()),'terminal_null_low':float(np.quantile(z[:,[0,19]].sum(1),.025)),'terminal_null_high':float(np.quantile(z[:,[0,19]].sum(1),.975))})
        summaries.extend({'scope':scope,'group':name,'proteins':n,'bin':j,'midpoint':(j+.5)/20,'observed_fraction':exp[j]+diff[j],'expected_fraction':exp[j],'difference':diff[j],'null_low':lo[j],'null_high':hi[j]} for j in range(20))
    sm=pd.DataFrame(summaries);tt=pd.DataFrame(tests);tt['q_BH_scope_family']=multipletests(tt.p_value,method='fdr_bh')[1]
    save(sm,out/f'{scope}_position_summary.tsv');save(tt,out/f'{scope}_group_omnibus.tsv');np.savez_compressed(out/f'{scope}_null_draws.npz',group_names=np.array(names),centered_protein_mean=np.array([null[i]/max(groups[n].sum(),1) for i,n in enumerate(names)]))
    return {'scope':scope,'proteins':len(d),'K':int(d.K.sum()),'Kla':int(d.Kla.sum()),'excluded':len(excluded),'draws':draws}

def main()->int:
    p=argparse.ArgumentParser();p.add_argument('--all-k',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--draws',type=int,default=999);p.add_argument('--seed',type=int,default=20260919);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    k=pd.read_csv(a.all_k,sep='\t',dtype={'taxon_id':str},keep_default_na=False);assert not k.duplicated(['protein_unit_id','k_position']).any()
    # P001 uses every primary case and censors known nonprimary Kla from the K denominator.
    prof=[]
    for pid,g0 in k.groupby('protein_unit_id',sort=True):
        g=g0.loc[~g0.is_censored.eq(1)];prof.append({'protein_unit_id':pid,'taxon_id':str(g.taxon_id.iloc[0]),'accession':g.accession.iloc[0],'n_K_evaluable':len(g),'n_primary_Kla':int(g.is_kla.sum()),'n_K_censored':int(g0.is_censored.sum()),'reported_Kla_fraction':float(g.is_kla.sum()/len(g))})
    save(pd.DataFrame(prof),a.output/'P001_PROTEIN_COVERAGE.tsv.gz')
    results=[analyze(k,s,*cols,a.output,a.draws,a.seed+i) for i,(s,cols) in enumerate(SCOPES.items())]
    (a.output/'P042_SUMMARY.json').write_text(json.dumps({'status':'PASS','results':results,'main_scope':'global_unknown','restricted_separate':True,'all_primary_completeness_mixed_sensitivity':True,'input_sha256':sha(a.all_k)},indent=2)+'\n');print(json.dumps(results,indent=2));return 0
