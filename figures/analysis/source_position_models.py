"""Between-source protein-relative position agreement under fixed K burdens."""
from __future__ import annotations
import argparse,json,hashlib,shutil,itertools
from pathlib import Path
import numpy as np,pandas as pd

def main()->None:
    """Compute a distinct CDF-position statistic, exact null means and selected null ranges."""
    p=argparse.ArgumentParser()
    for n in ['evidence','all-k','censored','reference-opportunity','components','output']:p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--draws',type=int,default=999);p.add_argument('--seed',type=int,default=20260919);a=p.parse_args();o=a.output;o.mkdir(parents=True,exist_ok=True);files=[a.evidence,a.all_k,a.censored,a.reference_opportunity,a.components,Path(__file__)];pd.DataFrame([dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files]).to_csv(o/'INPUT_MANIFEST.tsv',sep='\t',index=False)
    (o/'DESIGN.md').write_text('# P017 protein-relative position consistency\n\nExploratory, not preregistered. For every same-taxon unordered source-block pair, restrict to jointly reported proteins with exact complete reference lysine sets; censor current branch-specific known Kla positions, including selected/restricted/nonprimary-only and sequence-alias masks. Each source keeps its protein-specific reported burden. Divide full relative protein length into20 fixed equal bins; compare cumulative source fractions at the19 internal boundaries. Per-protein distance D=mean_j(Fa_j−Fb_j)², and pair distance is the equal-protein mean, not a pooled-site distance or P007 coordinate overlap. Under independent uniform subsets of a and b positions from N eligible K, E(D)=mean_j q_j(1−q_j)[(N−a)/(a(N−1))+(N−b)/(b(N−1))]; N=1 contributes0. This controls proteinK distribution and source burdens. All pairs have observed/expected distances;12 pairs at evenly spaced ranks of eligible shared protein count, tie by pairID, receive999 conditional-placement draws with a preselected seed. Their2.5–97.5% ranges are placement-null quantiles, not CI. No pairwise p/FDR or independent-experiment inference. Reuse component IDs accompany every pair. Positive E−D indicates closer positions than uniform placement; it does not establish biology versus detection causality.\n')
    e=pd.read_csv(a.evidence,sep='\t',usecols=['strict_site_id','block_id','taxon_id','primary_uniprot_accession','site_position'],dtype={'taxon_id':str});e=e[['strict_site_id','block_id','taxon_id','primary_uniprot_accession','site_position']].drop_duplicates();e['protein']=e.taxon_id+'|'+e.primary_uniprot_accession;k=pd.read_csv(a.all_k,sep='\t',dtype={'taxon_id':str},keep_default_na=False);c=set(pd.read_csv(a.censored,sep='\t').strict_site_id);op=pd.read_csv(a.reference_opportunity,sep='\t');valid=set(op.loc[op.full_reference_verified,'protein_unit_id']);k=k[k.protein_unit_id.isin(valid)&~k.strict_site_id.isin(c)];k['bin']=np.minimum(np.floor(20*(k.k_position-1)/(k.sequence_length-1)).astype(int),19)
    ku={pid:set(g.k_position) for pid,g in k.groupby('protein_unit_id')};kb={pid:np.bincount(g.bin,minlength=20) for pid,g in k.groupby('protein_unit_id')};posbin={pid:dict(zip(g.k_position,g.bin)) for pid,g in k.groupby('protein_unit_id')};raw={};data={};fail=[]
    for (block,pid),g in e.groupby(['block_id','protein']):
        sites=set(g.site_position);raw.setdefault(block,set()).add(pid)
        if pid not in ku or not sites<=ku[pid]:fail.append(dict(block_id=block,protein=pid,reported_sites=len(sites),reason='no_complete_unique_reference_or_nonK_position'));continue
        data[block,pid]=np.bincount([posbin[pid][x] for x in sites],minlength=20)
    pd.DataFrame(fail).to_csv(o/'SOURCE_PROTEIN_EXCLUSIONS.tsv',sep='\t',index=False)
    comp=pd.read_csv(a.components,sep='\t',dtype=str);comp=comp[comp.retained_in_primary.eq('True')&comp.scenario.eq('conservative_publication_reuse')];ev=pd.read_csv(a.evidence,sep='\t',usecols=['source_table_id','block_id']);bc=ev.merge(comp[['source_table_id','relationship_unit']],on='source_table_id').groupby('block_id').relationship_unit.agg(lambda x:sorted(set(x))[0]).to_dict();rows=[];contrib=[];pairprots={}
    for tax,g in e.groupby('taxon_id'):
        for ba,bb in itertools.combinations(sorted(g.block_id.unique()),2):
            common={x for x in raw[ba]&raw[bb] if x.startswith(tax+'|')};pp=sorted(x for x in common if (ba,x) in data and (bb,x) in data);pid=f'{tax}|{ba}|{bb}';observed=[];expected=[]
            for protein in pp:
                na=kb[protein];N=int(na.sum());aa=data[ba,protein];ab=data[bb,protein];n=int(aa.sum());m=int(ab.sum());fa=np.cumsum(aa)[:-1]/n;fb=np.cumsum(ab)[:-1]/m;q=np.cumsum(na)[:-1]/N;d=float(np.mean((fa-fb)**2));ex=float(np.mean(q*(1-q))*((N-n)/(n*(N-1))+(N-m)/(m*(N-1)))) if N>1 else 0.;observed.append(d);expected.append(ex);contrib.append(dict(pair_id=pid,taxon_id=tax,block_a=ba,block_b=bb,protein=protein,K=N,source_a_cases=n,source_b_cases=m,observed_CDF_distance=d,expected_CDF_distance=ex))
            r=dict(pair_id=pid,taxon_id=tax,block_a=ba,block_b=bb,component_a=bc[ba],component_b=bc[bb],same_known_component=bc[ba]==bc[bb],shared_reported_proteins=len(common),eligible_shared_proteins=len(pp),lost_common_proteins=len(common)-len(pp),status='estimable_descriptive' if pp else 'no_eligible_shared_protein')
            if pp:r.update(observed_distance=float(np.mean(observed)),expected_distance=float(np.mean(expected)),agreement_excess=float(np.mean(expected)-np.mean(observed)));pairprots[pid]=(ba,bb,pp)
            rows.append(r)
    pair=pd.DataFrame(rows);pair.to_csv(o/'PAIR_POSITION_CONSISTENCY.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(contrib).to_csv(o/'PROTEIN_POSITION_CONTRIBUTIONS.tsv.gz',sep='\t',index=False)
    selection=pair[pair.status.eq('estimable_descriptive')].sort_values(['eligible_shared_proteins','pair_id']);selected=selection.iloc[np.unique(np.linspace(0,len(selection)-1,12).astype(int))];rng=np.random.default_rng(a.seed);mc=[];draws=[]
    for r in selected.itertuples():
        ba,bb,pp=pairprots[r.pair_id];null=np.zeros(a.draws)
        for protein in pp:
            colors=kb[protein];aa=int(data[ba,protein].sum());bbn=int(data[bb,protein].sum());x=rng.multivariate_hypergeometric(colors,aa,size=a.draws);y=rng.multivariate_hypergeometric(colors,bbn,size=a.draws);df=np.cumsum(x,axis=1)[:,:-1]/aa-np.cumsum(y,axis=1)[:,:-1]/bbn;null+=np.mean(df**2,axis=1)/len(pp)
        se=null.std(ddof=1)/np.sqrt(a.draws);lo,hi=np.quantile(null,[.025,.975]);mc.append(dict(pair_id=r.pair_id,taxon_id=r.taxon_id,shared_proteins=len(pp),observed_distance=r.observed_distance,exact_expected_distance=r.expected_distance,MC_mean=null.mean(),MC_SE=se,mean_error_in_MCSE=(null.mean()-r.expected_distance)/se if se else 0.,null_q025=lo,null_q975=hi,draws=a.draws));draws.extend(dict(pair_id=r.pair_id,draw=j,null_distance=float(v)) for j,v in enumerate(null));print('MC',r.taxon_id,len(pp),flush=True)
    pd.DataFrame(mc).to_csv(o/'RULE_SELECTED_PAIR_NULL_RANGES.tsv',sep='\t',index=False);pd.DataFrame(draws).to_csv(o/'RULE_SELECTED_NULL_DRAWS.tsv.gz',sep='\t',index=False)
    summary=dict(pairs=len(pair),estimable=int(pair.status.eq('estimable_descriptive').sum()),protein_pair_rows=len(contrib),source_protein_exclusions=len(fail),MC_pairs=len(mc),draws=a.draws,max_MC_expectation_error_in_MCSE=max(abs(x['mean_error_in_MCSE']) for x in mc),same_component_pairs=int(pair.same_known_component.sum()));(o/'PAIR_SUMMARY.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
if __name__=='__main__':main()
