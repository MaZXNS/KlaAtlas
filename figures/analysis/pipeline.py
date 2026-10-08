"""Calculate the figure results from database inputs."""
from __future__ import annotations
import argparse, importlib, json, shutil, sqlite3, subprocess, sys, tempfile, time, platform
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

def invoke(module: str, args: list[str]) -> None:
    """Run a portable module's retained command line in the current interpreter."""
    old = sys.argv
    sys.argv = [module, *map(str, args)]
    try:
        importlib.import_module(module).main()
    finally:
        sys.argv = old

def save(d: pd.DataFrame, p: Path) -> None:
    """Write result tables."""
    p.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(p, sep='\t', index=False, na_rep='NA', compression={'method': 'gzip', 'mtime': 0} if p.suffix == '.gz' else None)

def basis_from_frozen(frame: pd.DataFrame, column: str, prefix: str, params: dict) -> tuple:
    """Recover the centered natural-spline transform from the stored design rows,
    knots and boundaries.
    """
    import patsy
    x = np.log1p(frame[column].to_numpy(float))
    params = dict(params)
    params['kn'] = tuple(params['kn'])
    design = patsy.dmatrix('cr(x, knots=kn, lower_bound=lo, upper_bound=hi) - 1', dict(x=x, **params))
    names = [c for c in frame if c.startswith(prefix)]
    full = np.asarray(design)
    target = frame[names].to_numpy(float)
    if np.linalg.matrix_rank(full) != full.shape[1]:
        raise ValueError('Frozen rows do not identify the original spline transform')
    transform = np.linalg.lstsq(full, target, rcond=None)[0]
    error = float(np.max(np.abs(full @ transform - target)))
    if error > 1e-10:
        raise ValueError(f'Cannot recover exact frozen centered spline: {error}')
    return (design.design_info, transform, names, params, error)

def peptide_analysis(inp: Path, out: Path) -> None:
    """Fit the three peptide-length models and the interval model."""
    import exact_conditional as exact, peptide_models as pm, patsy
    kernel = exact.ExactConditional(out / 'runtime')
    p = inp / 'inputs/position/peptide_localization'
    dest = out / 'peptide'
    dest.mkdir(parents=True, exist_ok=True)
    frames = {n: pd.read_csv(p / 'P032' / f'{n}_INFORMATIVE_ROWS.tsv.gz', sep='\t') for n in ['actual_shortest', 'actual_median', 'theoretical_shortest']}
    basis = basis_from_frozen(pd.concat(frames.values(), ignore_index=True), 'length_covariate', 'length_spline_', json.loads((p / 'P032/SPLINE_DESIGN.json').read_text()))
    design, T, names, params, error = basis
    curves = []
    for name, f in frames.items():
        model = pm.fitting(f, names + ['flank_KR_fraction', 'relative_position'], name, dest, exact, kernel)
        grid = np.linspace(f.length_covariate.min(), f.length_covariate.max(), 120)
        B = np.asarray(patsy.build_design_matrices([design], dict(x=np.log1p(np.r_[grid, 15]), **params))[0]) @ T
        D = B[:-1] - B[-1]
        for x, v in zip(grid, D):
            full = np.zeros(len(model['features']))
            for value, n in zip(v, names):
                full[model['features'].index(n)] = value
            be = full @ model['beta']
            se = np.sqrt(max(0, full @ model['cov'] @ full))
            curves.append(dict(model=name, x=x, reference=15, estimable=model['valid'], log2_relative_or=be / np.log(2), ci95_lower=(be - 1.96 * se) / np.log(2), ci95_upper=(be + 1.96 * se) / np.log(2)))
    save(pd.DataFrame(curves), dest / 'LENGTH_ASSOCIATION_CURVES.tsv')
    f = pd.read_csv(p / 'P034/within_peptide_INFORMATIVE_ROWS.tsv.gz', sep='\t')
    names = [c for c in f if c.startswith('end_spline_')]
    model = pm.fitting(f, names + ['signed_relative_peptide_position'], 'within_peptide', dest, exact, kernel)
    try:
        design, T, names, params, error = basis_from_frozen(f, 'symmetric_end_distance', 'end_spline_', json.loads((p / 'P034/SPLINE_DESIGN.json').read_text()))
        geom = []
        for length in [15, 25]:
            N = np.arange(length)
            C = length - 1 - N
            near = np.minimum(N, C)
            signed = (N - C) / (length - 1)
            B = np.asarray(patsy.build_design_matrices([design], dict(x=np.log1p(np.r_[near, (length - 1) / 2]), **params))[0]) @ T
            D = B[:-1] - B[-1]
            for i, v in enumerate(D):
                full = np.zeros(len(model['features']))
                for value, n in zip(v, names):
                    full[model['features'].index(n)] = value
                full[model['features'].index('signed_relative_peptide_position')] = signed[i]
                be = full @ model['beta']
                se = np.sqrt(max(0, full @ model['cov'] @ full))
                geom.append(dict(peptide_length=length, N_distance=N[i], C_distance=C[i], relative_N_position=N[i] / (length - 1), reference_N_position=(length - 1) / 2, log2_relative_or=be / np.log(2), ci95_lower=(be - 1.96 * se) / np.log(2), ci95_upper=(be + 1.96 * se) / np.log(2), display_role='geometry-consistent_model_profile_not_a_measured_peptide'))
        save(pd.DataFrame(geom), dest / 'GEOMETRY_CONSISTENT_END_PROFILES.tsv')
        (dest / 'SPLINE_RECOVERY.json').write_text(json.dumps({'algorithm': 'Centered spline basis transform', 'maximum_design_error': error, 'transform': T.tolist()}, indent=2))
    except Exception as err:
        raise
    obs = pd.read_csv(p / 'P035_ALL_TARGET_BOUND_SCORE_OBSERVATIONS.tsv.gz', sep='\t')
    units = []
    for key, g in obs.groupby(['source_table_id', 'score_type', 'peptide_key'], sort=True):
        r = g.iloc[0]
        values = np.array(sorted(set(g.score)))
        units.append(dict(source_table_id=key[0], score_type=key[1], peptide_key=key[2], peptide_K_count=r.peptide_K_count, nearest_other_K_distance=r.nearest_other_K_distance, all_observations_at_threshold=bool((values >= r.score_threshold).all())))
    u = pd.DataFrame(units)
    u['K_group'] = pd.cut(u.peptide_K_count, [0, 1, 2, 3, 5, np.inf], labels=['1', '2', '3', '4-5', '6+']).astype(str)
    u['gap_group'] = pd.cut(u.nearest_other_K_distance, [-1, 1, 3, 7, 15, np.inf], labels=['1', '2-3', '4-7', '8-15', '16+']).astype(str)
    u.loc[u.peptide_K_count.eq(1), 'gap_group'] = 'single_K_not_applicable'
    source = u.groupby(['source_table_id', 'score_type', 'K_group', 'gap_group'], observed=True).agg(peptide_units=('peptide_key', 'size'), threshold_pass=('all_observations_at_threshold', 'sum')).reset_index()
    source['threshold_fraction'] = source.threshold_pass / source.peptide_units
    pc = source[source.score_type.eq('probability_scalar')]
    display = pc.groupby(['K_group', 'gap_group'], observed=True).agg(source_tables=('source_table_id', 'nunique'), equal_table_threshold_fraction=('threshold_fraction', 'mean'), peptide_units=('peptide_units', 'sum')).reset_index()
    save(display, dest / 'P035_GEOMETRY_SCORE_CELLS.tsv')
    save(source, dest / 'P035_SOURCE_GEOMETRY_CELLS.tsv')

def sequence_analysis(inp: Path, out: Path) -> None:
    """Refit full matched sequence and paired symmetry families and exact histone models."""
    import ptm_logo as logo, sequence_contrasts as sc, sequence_symmetry as sym
    d = inp / 'inputs/figure_3'
    dest = out / 'sequence'
    dest.mkdir(parents=True, exist_ok=True)
    fg = pd.read_csv(d / 'a_foreground_a_CURRENT117983_SITE_WINDOWS.tsv.gz', sep='\t')
    bg = pd.read_csv(d / 'a_background_ALL_K_SEQUENCE_WINDOWS.tsv.gz', sep='\t')
    cnt, bcnt, h, stat = logo.build_ptm_logo(fg, bg)
    h.to_csv(dest / 'F4a_PTM_PREFERENCE_HEIGHTS.tsv', sep='\t')
    stat.to_csv(dest / 'F4a_PTM_POSITION_STATS.tsv', sep='\t', index=False)
    del fg, bg
    f = pd.read_csv(inp / 'inputs/sequence/shared/COMPLETE_WINDOW_FEATURES.tsv.gz', sep='\t', dtype={'taxon_id': str})
    f = f[f.global_is_censored.eq(0)].copy()
    f['is_kla'] = f.global_is_kla.astype(int)
    tax, meta, den = sc.sequence_effects(f, 'global_main_protein_position_decile')
    meta['estimable'] = True
    meta['reason'] = 'estimable'
    meta['family_size'] = 400
    meta['q_bh_400'] = multipletests(meta.p_value.fillna(1), method='fdr_bh')[1]
    save(meta, dest / 'b_global_main_P051_POSITION400.tsv')
    save(tax, dest / 'P051_TAXON.tsv.gz')
    f, *_ = sc.prepare(f)
    tax, summary, den = sym.run(f, 'global_main')
    save(summary, dest / 'c_d_global_main_P053_BH200.tsv')
    save(tax, dest / 'P053_TAXON.tsv.gz')
    invoke('histone_models', ['--opportunities', str(inp / 'inputs/histone/OPPORTUNITIES.tsv.gz'), '--output', str(dest / 'histone'), '--membership', 'recollected_supported'])

def functional_ac_analysis(inp: Path, out: Path) -> None:
    """Recount localization and refit stored compartment and RSA models."""
    import exact_conditional as exact, two_way_models as tw, structure_models as st
    from scipy.stats import norm
    dest = out / 'functional_ac'
    dest.mkdir(parents=True, exist_ok=True)
    d = inp / 'inputs/figure_4'
    ledger = pd.read_csv(inp / 'inputs/functional/HUMAN_SUBCELLULAR_PROTEIN_LEDGER.tsv.gz', sep='\t').fillna('')
    rows = []
    assigned = np.zeros(len(ledger), bool)
    with sqlite3.connect(f'file:{DATABASE}?mode=ro', uri=True) as conn:
        sites = pd.read_sql_query("SELECT accession,site_id FROM sites WHERE taxon_id='9606'", conn)
    sites['protein_unit_id'] = '9606|' + sites.accession.astype(str)
    counts = sites.groupby('protein_unit_id').site_id.nunique()
    for cat in ['Nucleus', 'Cytoplasm', 'Mitochondrion', 'ER', 'Golgi', 'Cell membrane', 'Secreted', 'Unknown']:
        mask = ~assigned if cat == 'Unknown' else ledger.categories.str.split(';').map(lambda x: cat in x).to_numpy()
        assigned |= mask
        rows.append(dict(category=cat, proteins=int(mask.sum()), sites_on_annotated_proteins=int(counts.reindex(ledger.loc[mask, 'protein_unit_id']).fillna(0).sum()), denominator_human_atlas_proteins=len(ledger), denominator_human_primary_sites=len(sites), fraction=float(mask.mean()), assignment='unassigned_no_qualifying_exact_reference_location' if cat == 'Unknown' else 'overlapping_multilabel_protein_annotation_not_site_localization_experiment'))
    save(pd.DataFrame(rows), dest / 'HUMAN_SUBCELLULAR_DISTRIBUTION.tsv')
    stack = pd.read_csv(inp / 'inputs/functional/models/p057/P057_CURRENT_SOURCE_CONDITIONAL_INPUT.tsv.gz', sep='\t')
    bindings = pd.read_csv(inp / 'inputs/functional/recovered_panel/shared/methods_batch5/annotation/EXACT_COMPARTMENT_BINDINGS.tsv', sep='\t')
    identity = stack[['original_protein_unit_id', 'taxon_id', 'sequence_sha256']].drop_duplicates()
    clusters = dict(zip(identity.original_protein_unit_id, identity.taxon_id.astype(str) + '|' + identity.sequence_sha256))
    kernel = exact.ExactConditional(dest / 'runtime')
    effects = []
    for cat in ['Nucleus', 'Cytoplasm', 'Mitochondrion', 'ER', 'Golgi', 'Cell membrane', 'Secreted']:
        ids = set(bindings.loc[bindings.category.eq(cat), 'protein_unit_id'])
        frame = stack[stack.original_protein_unit_id.isin(ids)].copy()
        meta, prepared, fit = tw.fit_two_way(frame, tw.BASE, kernel, clusters)
        for feature in tw.BASE[:4]:
            row = dict(category=cat, feature=feature, estimable=False, reason=meta['status'])
            if meta['status'] == 'estimable' and feature in meta['features']:
                j = meta['features'].index(feature)
                var = meta['covariance_two_way'][j][j]
                if var > 0:
                    b = meta['beta'][j] / np.log(2)
                    se = np.sqrt(var) / np.log(2)
                    row.update(estimable=True, reason='estimable', log2_or=b, standard_error=se, ci_low=b - 1.96 * se, ci_high=b + 1.96 * se, p_value=2 * norm.sf(abs(b / se)), source_components=meta['source_components'], protein_sequence_clusters=meta['sequence_protein_clusters'])
            effects.append(row)
    effect = pd.DataFrame(effects)
    effect['q_BH28'] = multipletests(effect.p_value.fillna(1), method='fdr_bh')[1]
    effect.loc[~effect.estimable, 'q_BH28'] = np.nan
    save(effect, dest / 'P057_COMPARTMENT_SEQUENCE_EFFECTS.tsv')
    f = pd.read_csv(d / 'P072_SOURCE_DATA.tsv.gz', sep='\t')
    absdiff = st.protein_balanced_difference(f, 1000, 20260922)
    save(absdiff, dest / 'P072_ABSOLUTE_DIFFERENCE.tsv')
    fitted90 = st.conditional_fit(f, ['residue_rsa'] + st.CONTROLS, 'global_unselected_pLDDT>=90')
    sensitivity = pd.read_csv(inp / 'inputs/functional/P072_PLDDT70_INPUT.tsv.gz', sep='\t')
    if sensitivity.duplicated(['protein_unit_id', 'site_position']).any():
        raise ValueError('Duplicated structure opportunity coordinate in pLDDT>=70 population')
    if not sensitivity.plddt.ge(70).all():
        raise ValueError('Unexpected pLDDT<70 in frozen sensitivity population')
    fitted70 = st.conditional_fit(sensitivity, ['residue_rsa'] + st.CONTROLS, 'global_unselected_pLDDT>=70')
    fitted = pd.concat([fitted90, fitted70], ignore_index=True)
    save(fitted, dest / 'P072_CONDITIONAL_EFFECTS.tsv')
    population = pd.DataFrame([dict(model=label, K=len(frame), Kla=int(frame.is_kla.sum()), proteins=frame.protein_unit_id.nunique()) for label, frame in [('global_unselected_pLDDT>=90', f), ('global_unselected_pLDDT>=70', sensitivity)]])
    save(population, dest / 'P072_POPULATIONS.tsv')

def resource_analysis(inp: Path, out: Path) -> None:
    """Derive stored atlas counts and exact 200-order source accumulation from core SQLite."""
    dest = out / 'resource'
    dest.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(f'file:{DATABASE}?mode=ro', uri=True) as conn:
        s = pd.read_sql_query('SELECT * FROM sites', conn)
        e = pd.read_sql_query('SELECT * FROM evidence', conn)
        p = pd.read_sql_query('SELECT * FROM studies', conn)
    e = e.merge(s[['site_id', 'taxon_id', 'accession', 'lysine_position', 'reference_sequence_id']], on='site_id', validate='many_to_one')
    summary = {}
    summary.update(sites=len(s), evidence=len(e), proteins=s[['taxon_id', 'accession']].drop_duplicates().shape[0], taxa=s.taxon_id.nunique(), publications=e.publication_id.nunique(), source_tables=e.source_table_id.nunique(), source_blocks=e.block_id.nunique())
    (dest / 'SUMMARY.json').write_text(json.dumps(summary, indent=2))
    tax = s.groupby('taxon_id').agg(species=('species', 'first'), sites=('site_id', 'size'), proteins=('accession', 'nunique')).join(e.groupby('taxon_id').agg(publications=('publication_id', 'nunique'), source_tables=('source_table_id', 'nunique'), source_blocks=('block_id', 'nunique'))).reset_index()
    save(tax, dest / 'TAXON_RESOURCE_COUNTS.tsv')
    yearmap = dict(zip(p.publication_id, p.year.astype(int)))
    ey = e[['site_id', 'publication_id']].drop_duplicates()
    ey['year'] = ey.publication_id.map(yearmap)
    first = ey.groupby('site_id').year.min()
    year = p.assign(year=p.year.astype(int)).groupby('year').size().rename('publications').to_frame().join(first.value_counts().rename('new_coordinates')).fillna(0).sort_index()
    year['cumulative_coordinates'] = year.new_coordinates.cumsum()
    save(year.reset_index(), dest / 'PUBLICATION_YEAR_AND_FIRST_COORDINATE.tsv')
    rng = np.random.default_rng(20260912)
    draws = []
    for taxon, ss in s.groupby('taxon_id', sort=True):
        ss = ss.sort_values('site_id')
        index = {x: i for i, x in enumerate(ss.site_id)}
        proteins, _ = pd.factorize(ss.accession, sort=True)
        eg = e[e.taxon_id.eq(taxon)]
        sets = {block: np.array(sorted({index[x] for x in g.site_id}), int) for block, g in eg.groupby('block_id')}
        blocks = sorted(sets)
        for it in range(200):
            seen = np.zeros(len(ss), bool)
            seenp = np.zeros(proteins.max() + 1, bool)
            newp = oldp = repeated = 0
            for step, j in enumerate(rng.permutation(len(blocks)), 1):
                ids = sets[blocks[j]]
                new = ~seen[ids]
                np_mask = ~seenp[proteins[ids]]
                newp += int(np.sum(new & np_mask))
                oldp += int(np.sum(new & ~np_mask))
                repeated += int(np.sum(~new))
                seen[ids] = True
                seenp[proteins[ids]] = True
                draws.append(dict(taxon_id=taxon, iteration=it, step=step, cumulative_new_protein_coordinates=newp, cumulative_existing_protein_coordinates=oldp, cumulative_added_source_support=repeated, union_sites=int(seen.sum())))
    raw = pd.DataFrame(draws)
    long = raw.melt(id_vars=['taxon_id', 'iteration', 'step'], var_name='measure')
    summary = long.groupby(['taxon_id', 'step', 'measure']).value.agg(median='median', mean='mean', low=lambda x: x.quantile(0.05), high=lambda x: x.quantile(0.95)).reset_index()
    save(summary, dest / 'SOURCE_ORDER_DECOMPOSITION.tsv')
    inputs = inp / 'inputs/resource/p017_inputs'
    all_ev = pd.read_csv(inp / 'inputs/resource/PRIMARY_EVIDENCE_BOUND.tsv.gz', sep='\t', low_memory=False)
    global_ev = all_ev[all_ev.design_class.eq('global_discovery') & all_ev.coverage_status.ne('selected_subset')]
    global_input = out / 'GLOBAL_MODELLING_EVIDENCE_BOUND.tsv.gz'
    save(global_ev, global_input)
    invoke('source_position_models', ['--evidence', str(global_input), '--all-k', str(inputs / 'ALL_K_REFERENCE.tsv.gz'), '--censored', str(inputs / 'CENSORED_GLOBAL_K.tsv'), '--reference-opportunity', str(inputs / 'REFERENCE_ELIGIBILITY.tsv'), '--components', str(inputs / 'COMPONENTS_CURRENT_CONSERVATIVE.tsv'), '--output', str(dest / 'p017')])
    global_input.unlink()
