"""Calculate all study- and namespace-specific enrichment models."""
from pathlib import Path
import pandas as pd
import numpy as np
import enrichment_models as model

def run(inp: Path, output: Path) -> dict:
    """Run the original 13-study models and preserve all eligible tests and exclusions."""
    output.mkdir(parents=True, exist_ok=True)
    specs = model.load_specs(inp / 'inputs/enrichment/populations', inp / 'inputs/figure_4/DISPLAY_208_UPDATED.tsv')
    populations = pd.concat([model.normalize_population(s) for s in specs], ignore_index=True)
    model.write_tsv(populations, output / 'ALL_POPULATION_ROWS.tsv.gz', gzip_output=True)
    model.write_tsv(populations[~populations.exclusion_reason.eq('eligible')], output / 'EXCLUDED_POPULATION_ROWS.tsv.gz', gzip_output=True)
    accession = pd.read_csv(inp / 'inputs/enrichment/memberships/ACCESSION_TERM_MEMBERSHIP.tsv.gz', sep='\t', dtype=str, keep_default_na=False)
    genes = pd.read_csv(inp / 'inputs/enrichment/memberships/GENE_TERM_MEMBERSHIP.tsv.gz', sep='\t', dtype=str, keep_default_na=False)
    catalog=pd.read_csv(inp/'inputs/enrichment/memberships/TERM_CATALOG.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
    results = []; eligibility = []; denominators = []; summary = []
    for spec in specs:
        pop = populations[populations.study.eq(spec.name)].copy()
        members = accession if spec.analysis_unit == 'accession' else genes
        r, e, d = model.analyse_study(pop, members, 10, 500, catalog)
        results.append(r); eligibility.append(e); denominators.append(d)
        mask = pop.exclusion_reason.eq('eligible')
        summary.append(dict(study=spec.name, analysis_unit=spec.analysis_unit,
                            frozen_population_rows=len(pop), covariate_complete_units=int(mask.sum()),
                            covariate_complete_foreground=int(pop.loc[mask, 'foreground_kla'].sum()),
                            excluded_units=int((~mask).sum())))
        print(f'{spec.name}: {len(r)} tests', flush=True)
    result = pd.concat(results, ignore_index=True)
    eligible = pd.concat(eligibility, ignore_index=True)
    den = pd.concat(denominators, ignore_index=True)
    model.write_tsv(result, output / 'ALL_TERM_RESULTS.tsv.gz', gzip_output=True)
    model.write_tsv(eligible, output / 'ALL_TERM_ELIGIBILITY.tsv.gz', gzip_output=True)
    model.write_tsv(den, output / 'NAMESPACE_DENOMINATORS.tsv')
    model.write_tsv(pd.DataFrame(summary), output / 'POPULATION_DENOMINATOR_SUMMARY.tsv')
    display = model.build_display(inp / 'inputs/figure_4/DISPLAY_208_UPDATED.tsv', result, eligible)
    model.write_tsv(display, output / 'DISPLAY_208_UPDATED.tsv')
    selected = display[['ontology', 'namespace', 'term_id', 'term_name']].drop_duplicates().copy()
    for ontology in ['GO', 'KEGG']:
        subset = selected[selected.ontology.eq(ontology)].copy()
        subset.insert(0, 'selected_order', np.arange(1, len(subset)+1))
        model.write_tsv(subset, output / f'SELECTED_{ontology}_TERMS.tsv')
    recurrence = []
    for (namespace, concept, name), g in display.groupby(['namespace', 'concept_id', 'term_name'], sort=False):
        values = pd.to_numeric(g.adjusted_or, errors='coerce'); q = pd.to_numeric(g.adjusted_q, errors='coerce')
        recurrence.append(dict(namespace=namespace, concept_id=concept, term_name=name,
                               estimable_studies=int(values.notna().sum()), positive_bh05=int(((values>1)&(q<.05)).sum()),
                               negative_bh05=int(((values<1)&(q<.05)).sum())))
    model.write_tsv(pd.DataFrame(recurrence), output / 'DISPLAY_RECURRENCE.tsv')
    counts = dict(studies=len(specs), bh_families=len(result[['study','namespace']].drop_duplicates()),
                  covariate_complete_units=int(populations.exclusion_reason.eq('eligible').sum()),
                  foreground_units=int(populations.loc[populations.exclusion_reason.eq('eligible'),'foreground_kla'].sum()),
                  all_tests=len(result), GO_tests=int(result.namespace.str.startswith('GO').sum()),
                  KEGG_tests=int(result.namespace.eq('KEGG').sum()),
                  GO_q_lt_005=int((result.namespace.str.startswith('GO')&result.adjusted_q.lt(.05)).sum()),
                  KEGG_q_lt_005=int((result.namespace.eq('KEGG')&result.adjusted_q.lt(.05)).sum()), display_cells=len(display))
    return counts
