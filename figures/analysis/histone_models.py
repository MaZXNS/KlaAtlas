"""KlaAtlas analysis routines."""
from __future__ import annotations

import argparse

import hashlib

import importlib.util

import json

import platform

import sys

from pathlib import Path

from typing import Any

import numpy as np

import pandas as pd

from scipy.stats import norm

from statsmodels.stats.multitest import multipletests

FEATURES=['kr_fraction','de_fraction','hydrophobicity','within_decile_position','theory_available','theory_length_log1p']

REPORTED=FEATURES[:4]

GROUPS=['verified_histone','annotated_other_protein']

def load(path: Path) -> Any:
    """Import the frozen original exact-conditional implementation."""
    spec=importlib.util.spec_from_file_location('histone_recollection_exact',path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    return module

def fit_group(frame: pd.DataFrame, group: str, features: list[str], exact: Any, kernel: Any) -> tuple[list[dict], dict, pd.DataFrame]:
    """Fit one group; retain original diagnostics and NA rather than unstable values."""
    data=exact.prepare_data(frame,features)
    fit=exact.fit_conditional(kernel,data,maxiter=300)
    used=frame.loc[data.row_index].copy()
    ok=bool(fit['success'] and fit['gradient_per_stratum_max']<1e-5 and fit['hessian_min_eigenvalue']>1e-7 and np.abs(fit['beta']).max()<20)
    rows=[]
    for feature in REPORTED:
        row=dict(group=group,feature=feature,estimable=False,reason='failed_diagnostic',log2_or=np.nan,standard_error=np.nan,ci_low=np.nan,ci_high=np.nan,p_value=np.nan)
        if ok and feature in data.features:
            i=data.features.index(feature);b=fit['beta'][i]/np.log(2);se=np.sqrt(fit['covariance_cluster'][i,i])/np.log(2)
            if np.isfinite(se) and se>0:
                row.update(estimable=True,reason='estimable',log2_or=b,standard_error=se,ci_low=b-1.96*se,ci_high=b+1.96*se,p_value=2*norm.sf(abs(b/se)))
        rows.append(row)
    crosstab=pd.crosstab(used.is_kla,used.theory_available_raw).reindex(index=[0,1],columns=[0,1],fill_value=0)
    meta=dict(group=group,status='estimable' if ok else 'failed_diagnostic',K=len(used),Kla=int(used.is_kla.sum()),sequence_units=len(data.protein_ids),strata=len(data.starts)-1,features=data.features,beta=fit['beta'].tolist(),covariance=fit['covariance_cluster'].tolist(),optimizer_success=fit['success'],optimizer_message=fit['message'],iterations=fit['iterations'],gradient_per_stratum=fit['gradient_per_stratum_max'],hessian_min=fit['hessian_min_eigenvalue'],theory_available_2x2=crosstab.to_dict())
    return rows,meta,used

def main() -> None:
    """Write all original and overlap results without selecting by significance."""
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--opportunities',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--exact-module',type=Path,default=Path(__file__).with_name('exact_conditional.py'))
    ap.add_argument('--membership',required=True)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(a.opportunities,sep='\t');assert not d.duplicated(['sequence_unit','k_position']).any()
    exact=load(a.exact_module);kernel=exact.ExactConditional(a.output/'runtime')
    results=[];diagnostics=[];flows=[]
    for scenario,features in [('original_six',FEATURES),('detectability_overlap', [x for x in FEATURES if x!='theory_available'])]:
        parent=d if scenario=='original_six' else d.loc[d.theory_available_raw.eq(1)].copy()
        entries=[]
        for group in GROUPS:
            frame=parent.loc[parent.histone_status.eq(group)].copy()
            rows,meta,used=fit_group(frame,group,features,exact,kernel)
            entries.extend(rows);meta.update(scenario=scenario,membership=a.membership);diagnostics.append(meta)
            used.to_csv(a.output/f'{scenario}__{group}_INFORMATIVE.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
            flows.append(dict(membership=a.membership,scenario=scenario,group=group,parent_K=len(d[d.histone_status.eq(group)]),eligible_K=len(frame),eligible_cases=int(frame.is_kla.sum()),eligible_units=frame.sequence_unit.nunique(),informative_K=meta['K'],informative_cases=meta['Kla'],informative_units=meta['sequence_units'],informative_strata=meta['strata'],status=meta['status']))
            print(a.membership,scenario,group,meta['status'],meta['sequence_units'],meta['Kla'],flush=True)
        estimates=pd.DataFrame(entries)
        for feature in REPORTED:
            v=estimates[estimates.feature.eq(feature)].set_index('group')
            row=dict(group='histone_minus_other',feature=feature,estimable=False,reason='group_nonestimable',log2_or=np.nan,standard_error=np.nan,ci_low=np.nan,ci_high=np.nan,p_value=np.nan)
            if v.estimable.all():
                b=v.loc[GROUPS[0],'log2_or']-v.loc[GROUPS[1],'log2_or'];se=np.sqrt((v.standard_error**2).sum())
                row.update(estimable=True,reason='estimable',log2_or=b,standard_error=se,ci_low=b-1.96*se,ci_high=b+1.96*se,p_value=2*norm.sf(abs(b/se)))
            entries.append(row)
        table=pd.DataFrame(entries);assert len(table)==12
        table['q_BH12']=multipletests(table.p_value.fillna(1),method='fdr_bh')[1];table.loc[~table.estimable,'q_BH12']=np.nan
        table['membership']=a.membership;table['scenario']=scenario
        table['annotation_evidence_scope']=('reviewed_unreviewed_and_historical_exact_supported' if a.membership=='recollected_archived_supported' else 'reviewed_and_unreviewed_exact_supported') if a.membership.endswith('supported') else 'reviewed_supported'
        table['group_display']=table.group.map({'verified_histone':'Histone group','annotated_other_protein':'Other annotated proteins','histone_minus_other':'Histone minus other'})
        results.append(table)
    pd.concat(results,ignore_index=True).to_csv(a.output/'EFFECTS.tsv',sep='\t',index=False)
    pd.DataFrame(flows).to_csv(a.output/'DENOMINATORS.tsv',sep='\t',index=False)
    (a.output/'DIAGNOSTICS.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    meta={'input':str(a.opportunities),'input_sha256':hashlib.sha256(a.opportunities.read_bytes()).hexdigest(),'exact_module':str(a.exact_module),'exact_module_sha256':hashlib.sha256(a.exact_module.read_bytes()).hexdigest(),'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'membership':a.membership,'fixed_BH_family_per_scenario':12,'scaling':'unchanged frozen parent','contrast_SE':'sqrt(var_histone+var_other), disjoint groups'}
    (a.output/'RUN.json').write_text(json.dumps(meta,indent=2)+'\n')
