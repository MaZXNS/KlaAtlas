#!/usr/bin/env python3
"""Run KlaAtlas analyses and draw figures from the single frozen database."""
from __future__ import annotations
import argparse,importlib,json,shutil,sqlite3,sys,tempfile,subprocess,time
from datetime import datetime,timezone
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/'analysis'));sys.path.insert(0,str(ROOT/'code'))
from database_inputs import restore
import pipeline


def run_analysis(database:Path,inputs:Path,output:Path,figures:list[str]) -> None:
    """Calculate each selected figure's results using its original frozen populations."""
    pipeline.DATABASE=database
    if '1' in figures:pipeline.resource_analysis(inputs,output)
    if '2' in figures:
        pipeline.peptide_analysis(inputs,output)
        pipeline.invoke('position_models',['--all-k',str(inputs/'inputs/sequence/shared/ALL_K_CURRENT_PRIMARY.tsv.gz'),'--output',str(output/'position')])
        pipeline.invoke('clustering_models',['--all-k',str(inputs/'inputs/sequence/shared/ALL_K_CURRENT_PRIMARY.tsv.gz'),'--output',str(output/'clustering')])
    if '3' in figures:pipeline.sequence_analysis(inputs,output)
    if '4' in figures:
        pipeline.functional_ac_analysis(inputs,output)
        import enrichment
        counts=enrichment.run(inputs,output/'functional_de');print(json.dumps(counts),flush=True)


def table(source:Path,destination:Path) -> None:
    """Place one required result in a temporary drawing directory; never silently substitute."""
    if not source.is_file():raise FileNotFoundError(f'Required analysis result not produced: {source.name}')
    destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,destination)


def drawing_inputs(inputs:Path,analysis:Path|None,destination:Path,figures:list[str]) -> None:
    """Use generated results for all scientific panels, or explicitly select database results."""
    mapping={
      '1':[('counts.json','inputs/figure_1/01_SUMMARY.json','resource/SUMMARY.json'),('taxon_counts.tsv','inputs/figure_1/02_TAXON_RESOURCE_COUNTS.tsv','resource/TAXON_RESOURCE_COUNTS.tsv'),('publication_years.tsv','inputs/figure_1/03_PUBLICATION_YEAR_AND_FIRST_COORDINATE.tsv','resource/PUBLICATION_YEAR_AND_FIRST_COORDINATE.tsv'),('source_accumulation.tsv','inputs/figure_1/04_SOURCE_ORDER_DECOMPOSITION.tsv','resource/SOURCE_ORDER_DECOMPOSITION.tsv'),('position_agreement.tsv','inputs/figure_1/05_PAIR_POSITION_CONSISTENCY.tsv','resource/p017/PAIR_POSITION_CONSISTENCY.tsv')],
      '2':[('coverage.tsv.gz','inputs/figure_2/01_P001_PROTEIN_COVERAGE.tsv.gz','position/P001_PROTEIN_COVERAGE.tsv.gz'),('positions.tsv','inputs/figure_2/02_global_unknown_position_summary.tsv','position/global_unknown_position_summary.tsv'),('terminal_tests.tsv','inputs/figure_2/03_global_unknown_group_omnibus.tsv','position/global_unknown_group_omnibus.tsv'),('clustering.tsv','inputs/figure_2/04_P045_BURDEN_MULTISCALE_BOOTSTRAP.tsv','clustering/P045_BURDEN_MULTISCALE_BOOTSTRAP.tsv'),('length_profiles.tsv','inputs/figure_2/05_LENGTH_ASSOCIATION_CURVES.tsv','peptide/LENGTH_ASSOCIATION_CURVES.tsv'),('peptide_geometry.tsv','inputs/figure_2/06_GEOMETRY_CONSISTENT_END_PROFILES.tsv','peptide/GEOMETRY_CONSISTENT_END_PROFILES.tsv'),('localization.tsv','inputs/figure_2/07_P035_GEOMETRY_SCORE_CELLS.tsv','peptide/P035_GEOMETRY_SCORE_CELLS.tsv')],
      '3':[('logo.tsv','inputs/figure_3/F4a_PTM_PREFERENCE_HEIGHTS.tsv','sequence/F4a_PTM_PREFERENCE_HEIGHTS.tsv'),('sequence_effects.tsv','inputs/figure_3/b_global_main_P051_POSITION400.tsv','sequence/b_global_main_P051_POSITION400.tsv'),('symmetry.tsv','inputs/figure_3/c_d_global_main_P053_BH200.tsv','sequence/c_d_global_main_P053_BH200.tsv'),('histone_effects.tsv','inputs/figure_3/e_recollected_supported_EFFECTS.tsv','sequence/histone/EFFECTS.tsv'),('histone_denominators.tsv','inputs/figure_3/e_recollected_supported_DENOMINATORS.tsv','sequence/histone/DENOMINATORS.tsv')],
      '4':[('localizations.tsv','inputs/figure_4/SourceData_4a_subcellular.tsv','functional_ac/HUMAN_SUBCELLULAR_DISTRIBUTION.tsv'),('compartment_effects.tsv','inputs/figure_4/SourceData_4b_P057.tsv','functional_ac/P057_COMPARTMENT_SEQUENCE_EFFECTS.tsv'),('rsa_difference.tsv','inputs/figure_4/P072_ABSOLUTE_DIFFERENCE.tsv','functional_ac/P072_ABSOLUTE_DIFFERENCE.tsv'),('rsa_effects.tsv','inputs/figure_4/P072_CONDITIONAL_EFFECTS.tsv','functional_ac/P072_CONDITIONAL_EFFECTS.tsv'),('enrichment.tsv','inputs/figure_4/DISPLAY_208_UPDATED.tsv','functional_de/DISPLAY_208_UPDATED.tsv')]
    }
    for figure in figures:
        target=destination/('F'+figure);target.mkdir(parents=True,exist_ok=True)
        for name,frozen,calculated in mapping.get(figure,[]):table(analysis/calculated if analysis is not None else inputs/frozen,target/name)
        if figure=='1':table(inputs/'inputs/figure_1/Figure_1a_SOURCE_TEXT.json',target/'selection.json')
        if figure=='4':
            # RSA values are the same frozen model covariates used by the analysis;
            # binning is a drawing step, not a substitute for the statistical fit.
            frame=pd.read_csv(inputs/'inputs/figure_4/P072_SOURCE_DATA.tsv.gz',sep='\t',usecols=['residue_rsa','is_kla']);edges=np.linspace(0,1.25,36);rows=[]
            for label,group in frame.groupby('is_kla'):
                density,_=np.histogram(group.residue_rsa,edges,density=True)
                rows.extend(dict(reported_Kla=label,bin_lower=edges[i],bin_upper=edges[i+1],density=value) for i,value in enumerate(density))
            pd.DataFrame(rows).to_csv(target/'rsa_histogram.tsv',sep='\t',index=False)
            cells=pd.read_csv(target/'enrichment.tsv',sep='\t',dtype={'term_id':str})
            for ontology,name in [('GO','go_terms.tsv'),('KEGG','kegg_terms.tsv')]:cells.loc[cells.ontology.eq(ontology),['ontology','namespace','term_id','term_name']].drop_duplicates().to_csv(target/name,sep='\t',index=False)
        if figure=='5':
            path=inputs/'inputs/figure_5/live_public_20260922';tracks=json.loads((path/'plot_tracks.json').read_text());meta=json.loads((path/'api_meta.json').read_text());meta=meta.get('data',meta)
            (target/'tracks.json').write_text(json.dumps(tracks));(target/'counts.json').write_text(json.dumps(meta['counts']))


def draw(figure:str,data:Path,output:Path) -> None:
    """Draw the existing layouts, retaining fixed workflow and interface material."""
    output.mkdir(parents=True,exist_ok=True)
    if figure=='1':
        with tempfile.TemporaryDirectory(prefix='klaatlas_panels_') as directory:
            import figure1a,figure1bf,pymupdf
            work=Path(directory);figure1a.draw(data/'F1',work);figure1bf.render_panels_bf(data/'F1',work);combined=pymupdf.open()
            for stem in ['Figure_1a','Figure_1b_f']:
                with pymupdf.open(work/(stem+'.pdf')) as original:
                    if original[0].rect.width>504.001:
                        factor=504/original[0].rect.width;page=combined.new_page(width=504,height=original[0].rect.height*factor);page.show_pdf_page(page.rect,original,0)
                    else:combined.insert_pdf(original)
            combined.save(output/'Figure_1.pdf',garbage=4,deflate=True);combined.close()
    elif figure=='S1':
        source=ROOT/'images/Figure_S1.pdf';target=output/'Figure_S1.pdf'
        if source.exists():
            if source.resolve()!=target.resolve():shutil.copy2(source,target)
        else:
            import pymupdf
            with pymupdf.open(ROOT/'assets/Figure_S1.svg') as document:target.write_bytes(document.convert_to_pdf())
    else:importlib.import_module('figure'+figure).draw(data/('F'+figure),output)


def main() -> None:
    """Run analysis, plot frozen results explicitly, or execute both stages."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['analysis','plot','all'],default='all')
    parser.add_argument('--figure',choices=['all','1','2','3','4','5','S1'],default='all')
    parser.add_argument('--database',type=Path,default=ROOT.parent/'database/klaatlas.sqlite')
    parser.add_argument('--output',type=Path,default=Path('results'))
    parser.add_argument('--analysis-dir',type=Path,help='With stage plot, use existing calculated results instead of database results')
    args=parser.parse_args();database=args.database.resolve();output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    if args.analysis_dir is not None and args.stage!='plot':parser.error('--analysis-dir is used only with --stage plot')
    figures=['1','2','3','4','5','S1'] if args.figure=='all' else [args.figure]
    scientific_figures=[figure for figure in figures if figure not in {'5','S1'}]
    started=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='klaatlas_inputs_') as directory:
        inputs=Path(directory);manifest=restore(database,inputs)
        if args.stage!='plot':
            if scientific_figures:run_analysis(database,inputs,output/'analysis',figures)
            else:print('Selected figures are fixed visual material and have no statistical analysis stage.',flush=True)
        if args.stage!='analysis':
            if not scientific_figures:print('Drawing fixed visual material; no statistical analysis.',flush=True)
            elif args.stage=='plot':print('Drawing explicitly from existing calculated results; no model fitting.' if args.analysis_dir else 'Drawing explicitly from frozen database results; no model fitting.',flush=True)
            else:print('Drawing from the results generated by this run; no result fallback.',flush=True)
            data=inputs/'drawing';drawing_inputs(inputs,output/'analysis' if args.stage=='all' else args.analysis_dir.resolve() if args.analysis_dir else None,data,figures)
            for figure in figures:draw(figure,data,output/'figures');print('Figure '+figure+' complete',flush=True)
    # Scientific provenance accompanies derived results, without a separate test tool.
    import importlib.metadata
    parameters={'F1':{'source_orders':200,'seed':20260912,'pair_null_draws':999,'pair_seed':20260919},'F2':{'position_draws':999,'position_seed':20260919,'clustering_draws_per_protein':10000,'clustering_seed':20260904,'bootstrap_replicates':1000,'bootstrap_seed':20260919,'window_widths':[10,20,30,50,100]},'F3':{'sequence_BH_family':400,'paired_contrast_BH_family':200,'histone_BH_family_per_scenario':12},'F4':{'compartment_BH_family':28,'RSA_bootstrap':1000,'RSA_bootstrap_seed':20260922,'RSA_conditional_method':'R survival::clogit Efron','RSA_thresholds':[90,70],'enrichment_term_size':[10,500],'enrichment_BH':'within study x annotation namespace'}}
    environment={'python':sys.version.split()[0],**{name:importlib.metadata.version(name) for name in ['numpy','pandas','scipy','statsmodels','matplotlib','patsy','logomaker','PyMuPDF']}}
    (output/'run_metadata.json').write_text(json.dumps({'stage':args.stage,'figures':figures,'scientific_result_source':'fixed_material' if not scientific_figures else 'generated_analysis' if args.stage=='all' else 'previous_generated_analysis' if args.analysis_dir else 'frozen_database_results' if args.stage=='plot' else 'generated_analysis','fixed_material':['selection workflow','S1 workflow','website interface'],'input_versions':manifest,'parameters':parameters,'environment':environment,'completed_UTC':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-started},indent=2))

if __name__=='__main__':main()
