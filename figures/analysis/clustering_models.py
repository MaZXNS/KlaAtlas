"""KlaAtlas analysis routines."""
from __future__ import annotations

import argparse, hashlib, importlib.util, json, sys

from pathlib import Path

import numpy as np

import pandas as pd

from scipy.stats import trim_mean

from statsmodels.stats.multitest import multipletests

SCOPES={"global_unknown":("global_is_kla","global_is_censored")}

WIDTHS=[10,20,30,50,100]

METRICS=['centered_log2_nearest']+[f'fraction_excess_{w}' for w in WIDTHS]


def save(d:pd.DataFrame,p:Path)->None:d.to_csv(p,sep='\t',index=False,compression={'method':'gzip','mtime':0} if p.suffix=='.gz' else None)

def main()->int:
    p=argparse.ArgumentParser();p.add_argument('--all-k',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--draws',type=int,default=10000);p.add_argument('--bootstrap',type=int,default=1000);p.add_argument('--seed',type=int,default=20260904);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    import multiscale_kernel as v3
    kernel=v3.MultiscaleKernel(a.output/'runtime');k=pd.read_csv(a.all_k,sep='\t',dtype={'taxon_id':str},keep_default_na=False);assert not k.duplicated(['protein_unit_id','k_position']).any();summ=[]
    for si,(scope,(case_col,censor_col)) in enumerate(SCOPES.items()):
        dest=a.output/scope;dest.mkdir(exist_ok=True);rows=[];excluded=[];nullsum=np.zeros((a.draws,len(METRICS)))
        for i,(pid,g0) in enumerate(k.groupby('protein_unit_id',sort=True)):
            g=g0.loc[~g0[censor_col].eq(1)];pos=g.k_position.to_numpy(int);cases=g.loc[g[case_col].eq(1),'k_position'].to_numpy(int)
            if len(cases)<2 or len(cases)>=len(pos):excluded.append({'protein_unit_id':pid,'K':len(pos),'Kla':len(cases),'reason':'fewer_than_two_cases' if len(cases)<2 else 'no_comparator_K'});continue
            seed=(a.seed+si*1000003+int(hashlib.sha256(pid.encode()).hexdigest()[:12],16))%(2**32);vals,draw=kernel.evaluate(pos,cases,WIDTHS,a.draws,seed);nullsum+=draw;rows.append({'protein_unit_id':pid,'taxon_id':str(g.taxon_id.iloc[0]),'accession':g.accession.iloc[0],'k_count':len(pos),'kla_count':len(cases),**vals})
            if i and i%1000==0:
                save(pd.DataFrame(rows),dest/'protein_statistics.checkpoint.tsv.gz');np.save(dest/'global_null_sum.checkpoint.npy',nullsum);(dest/'CHECKPOINT.json').write_text(json.dumps({'processed_protein_groups':i,'estimable':len(rows),'draws':a.draws})+'\n');print(scope,'P045 groups',i,'estimable',len(rows),flush=True)
        d=pd.DataFrame(rows);N=len(d);q=multipletests(np.r_[d.p_nearest,d.p_scan_family],method='fdr_bh')[1];d['q_nearest_joint_bh']=q[:N];d['q_scan_joint_bh']=q[N:];save(d,dest/'protein_statistics.tsv.gz');save(pd.DataFrame(excluded),dest/'excluded_proteins.tsv.gz');null=nullsum/N;save(pd.DataFrame(null,columns=METRICS),dest/'global_null.tsv.gz')
        gt=[]
        for j,m in enumerate(METRICS):
            eff=float(d[m].mean());pv=float((1+((null[:,j]<=eff) if j==0 else (null[:,j]>=eff)).sum())/(a.draws+1));gt.append({'metric':m,'effect':eff,'p_value':pv})
        gt=pd.DataFrame(gt);gt['q_BH6']=multipletests(gt.p_value,method='fdr_bh')[1];save(gt,dest/'global_tests.tsv')
        if scope=='global_unknown':
            rng=np.random.default_rng(20260919);bur=pd.cut(d.kla_count,[1,4,9,19,np.inf],labels=['2-4','5-9','10-19','20+']);br=[]
            for label in bur.cat.categories:
                g=d[bur.eq(label)]
                for w in WIDTHS:
                    x=g[f'fraction_excess_{w}'].to_numpy();boots=np.array([np.mean(x[ix]) for ix in rng.integers(0,len(x),(a.bootstrap,len(x)))])
                    br.append({'burden':str(label),'width':w,'proteins':len(x),'estimate':float(x.mean()),'ci95_low':float(np.quantile(boots,.025)),'ci95_high':float(np.quantile(boots,.975))})
            save(pd.DataFrame(br),a.output/'P045_BURDEN_MULTISCALE_BOOTSTRAP.tsv')
        summ.append({'scope':scope,'tested_proteins':N,'excluded':len(excluded),'Kla':int(d.kla_count.sum()),'draws':a.draws})
    (a.output/'P045_SUMMARY.json').write_text(json.dumps({'status':'PASS','results':summ,'bootstrap':a.bootstrap,'fixed_joint_family':'2 per protein across all tested proteins'},indent=2)+'\n');print(json.dumps(summ,indent=2));return 0
