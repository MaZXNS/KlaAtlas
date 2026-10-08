"""Frozen-population functional models with study-by-namespace BH correction."""

from __future__ import annotations

import json

import math

import warnings

from dataclasses import dataclass

from pathlib import Path

import numpy as np

import pandas as pd

import scipy.stats as st

import statsmodels.api as sm

from statsmodels.stats.multitest import multipletests

NAMESPACES = ("GO_BP", "GO_CC", "GO_MF", "KEGG")

GENE_STUDIES = {"Hypertrophic scar", "MM.1S WT–LenR"}

@dataclass(frozen=True)
class StudySpec:
    name: str
    species: str
    analysis_unit: str
    population_path: Path

def write_tsv(frame: pd.DataFrame, path: Path, *, gzip_output: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if gzip_output:
        frame.to_csv(path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    else:
        frame.to_csv(path, sep="\t", index=False)

def safe_name(value: str) -> str:
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in value)

def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})

def load_specs(source_run: Path, display_template: Path) -> list[StudySpec]:
    display = pd.read_csv(display_template, sep="\t", dtype=str, keep_default_na=False)
    ordered = display[["study", "species"]].drop_duplicates().itertuples(index=False)
    specs: list[StudySpec] = []
    for study, species in ordered:
        analysis_unit = "gene_id" if study in GENE_STUDIES else "accession"
        path = source_run / safe_name(study) / "population.tsv.gz"
        if not path.exists():
            raise FileNotFoundError(path)
        specs.append(StudySpec(study, species, analysis_unit, path))
    if len(specs) != 13:
        raise ValueError(f"expected 13 frozen studies, found {len(specs)}")
    return specs

def normalize_population(spec: StudySpec) -> pd.DataFrame:
    raw = pd.read_csv(spec.population_path, sep="\t", compression="gzip", dtype=str, keep_default_na=False)
    if spec.analysis_unit == "accession":
        out = pd.DataFrame(
            {
                "study": spec.name,
                "species": spec.species,
                "analysis_unit": "accession",
                "unit_id": raw["accession"].astype(str),
                "foreground_kla": as_bool(raw["foreground_kla"]),
                "abundance": pd.to_numeric(raw["abundance"], errors="coerce"),
                "k_count": pd.to_numeric(raw["k_count"], errors="coerce"),
                "source_accessions": raw["accession"].astype(str),
                "source_status": "frozen_unique_accession_row",
            }
        )
        out["exclusion_reason"] = np.where(
            out["abundance"].isna(),
            "missing_positive_abundance",
            np.where(out["k_count"].isna(), "missing_unambiguous_reference_K_count", "eligible"),
        )
    else:
        out = pd.DataFrame(
            {
                "study": spec.name,
                "species": spec.species,
                "analysis_unit": "gene_id",
                "unit_id": raw["gene_id"].astype(str),
                "foreground_kla": as_bool(raw["foreground_kla"]),
                "abundance": pd.to_numeric(raw["abundance_gene_sum"], errors="coerce"),
                "k_count": pd.to_numeric(raw["k_opportunity_max"], errors="coerce"),
                "source_accessions": raw["accessions"].astype(str),
                "source_status": raw["common_covariate_status"].astype(str),
            }
        )
        eligible = raw["common_covariate_status"].eq("eligible_unique_detected_accession")
        out["exclusion_reason"] = np.where(
            eligible,
            "eligible",
            "gene_unit_multiple_accessions_or_missing_abundance_or_reference_K",
        )
    if out["unit_id"].duplicated().any():
        raise ValueError(f"duplicate analysis units in {spec.name}")
    return out

def fisher_exact(a: int, b: int, c: int, d: int) -> dict[str, object]:
    table = np.asarray([[a, b], [c, d]], dtype=int)
    p_value = float(st.fisher_exact(table, alternative="two-sided").pvalue)
    try:
        estimate = st.contingency.odds_ratio(table, kind="conditional")
        if np.any(table == 0):
            return {
                "fisher_or": float(estimate.statistic),
                "fisher_ci_low": np.nan,
                "fisher_ci_high": np.nan,
                "fisher_p": p_value,
                "fisher_status": "conditional_exact_ci_unavailable_zero_cell",
            }
        interval = estimate.confidence_interval(0.95)
        return {
            "fisher_or": float(estimate.statistic),
            "fisher_ci_low": float(interval.low),
            "fisher_ci_high": float(interval.high),
            "fisher_p": p_value,
            "fisher_status": "conditional_exact",
        }
    except Exception as exc:  # pragma: no cover - defensive carrier
        return {
            "fisher_or": np.nan,
            "fisher_ci_low": np.nan,
            "fisher_ci_high": np.nan,
            "fisher_p": p_value,
            "fisher_status": f"failed:{type(exc).__name__}",
        }

def logistic(frame: pd.DataFrame) -> dict[str, object]:
    base = {
        "adjusted_or": np.nan,
        "adjusted_ci_low": np.nan,
        "adjusted_ci_high": np.nan,
        "adjusted_p": np.nan,
        "adjusted_status": "not_estimable",
    }
    if frame["term"].nunique() < 2 or frame["foreground_kla"].nunique() < 2:
        base["adjusted_status"] = "no_term_or_outcome_variation"
        return base
    two_by_two = pd.crosstab(frame["foreground_kla"], frame["term"]).reindex(
        index=[False, True], columns=[False, True], fill_value=0
    )
    if (two_by_two.to_numpy() == 0).any():
        base["adjusted_status"] = "complete_separation_in_term_2x2"
        return base
    design = frame[["term", "z_log1p_abundance", "z_log1p_k_count"]].astype(float)
    design = sm.add_constant(design, has_constant="add")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            fit = sm.Logit(frame["foreground_kla"].astype(int), design).fit(disp=False, maxiter=200)
        if not fit.mle_retvals.get("converged", False):
            base["adjusted_status"] = "nonconverged"
            return base
        beta = float(fit.params["term"])
        se = float(fit.bse["term"])
        if not np.isfinite(beta + se) or abs(beta) >= 20:
            base["adjusted_status"] = "separation_or_unstable"
            return base
        base.update(
            adjusted_or=math.exp(beta),
            adjusted_ci_low=math.exp(beta - 1.96 * se),
            adjusted_ci_high=math.exp(beta + 1.96 * se),
            adjusted_p=float(fit.pvalues["term"]),
            adjusted_status="estimable",
        )
    except Exception as exc:
        base["adjusted_status"] = f"fit_failed:{type(exc).__name__}"
    return base

def standardize(frame: pd.DataFrame, column: str) -> pd.Series:
    transformed = np.log1p(frame[column].astype(float))
    sd = float(transformed.std(ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError(f"constant or invalid {column}")
    return (transformed - transformed.mean()) / sd

def analyse_study(
    population_all: pd.DataFrame,
    membership: pd.DataFrame,
    min_term: int,
    max_term: int,
    term_catalog: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    study = str(population_all["study"].iloc[0])
    unit_type = str(population_all["analysis_unit"].iloc[0])
    eligible = population_all.loc[population_all["exclusion_reason"].eq("eligible")].copy()
    terms_all: list[dict[str, object]] = []
    results: list[dict[str, object]] = []
    namespace_rows: list[dict[str, object]] = []
    study_membership = membership.loc[membership["unit_id"].isin(set(eligible["unit_id"]))].copy()
    for namespace in NAMESPACES:
        ns_membership = study_membership.loc[study_membership["namespace"].eq(namespace)].drop_duplicates(
            ["unit_id", "term_id"]
        )
        annotated_ids = set(ns_membership["unit_id"])
        background = eligible.loc[eligible["unit_id"].isin(annotated_ids)].copy()
        if len(background):
            background["z_log1p_abundance"] = standardize(background, "abundance")
            background["z_log1p_k_count"] = standardize(background, "k_count")
        namespace_rows.append(
            {
                "study": study,
                "analysis_unit": unit_type,
                "namespace": namespace,
                "eligible_covariate_units": len(eligible),
                "annotated_model_units": len(background),
                "annotated_foreground": int(background["foreground_kla"].sum()) if len(background) else 0,
                "annotated_controls": int((~background["foreground_kla"]).sum()) if len(background) else 0,
            }
        )
        groups = {(str(term_id),str(term_name)):group for (term_id,term_name),group in ns_membership.groupby(["term_id","term_name"],sort=True)}
        if term_catalog is None:
            term_keys = sorted(groups)
        else:
            catalog = term_catalog.loc[term_catalog["analysis_unit"].eq(unit_type) & term_catalog["namespace"].eq(namespace),["term_id","term_name"]].drop_duplicates()
            term_keys = sorted((str(row.term_id),str(row.term_name)) for row in catalog.itertuples(index=False))
            missing=set(groups)-set(term_keys)
            if missing:raise ValueError("Mapped annotation terms are absent from the frozen term catalog")
        for term_id,term_name in term_keys:
            # The original eligibility inventory enumerates terms actually mapped
            # in this namespace population, not unobserved zero-overlap terms.
            if (term_id,term_name) not in groups:continue
            group=groups[term_id,term_name]
            members = set(group["unit_id"])
            term_size = int(background["unit_id"].isin(members).sum()) if len(background) else 0
            term_status = "eligible" if min_term <= term_size <= max_term else "outside_term_size_10_500"
            terms_all.append(
                {
                    "study": study,
                    "analysis_unit": unit_type,
                    "namespace": namespace,
                    "term_id": term_id,
                    "term_name": term_name,
                    "annotated_model_units": len(background),
                    "term_background_units": term_size,
                    "eligibility_status": term_status,
                }
            )
            if term_status != "eligible":
                continue
            frame = background.copy()
            frame["term"] = frame["unit_id"].isin(members)
            a = int((frame["foreground_kla"] & frame["term"]).sum())
            b = int((frame["foreground_kla"] & ~frame["term"]).sum())
            c = int((~frame["foreground_kla"] & frame["term"]).sum())
            d = int((~frame["foreground_kla"] & ~frame["term"]).sum())
            foreground_total = a + b
            control_total = c + d
            row = {
                "study": study,
                "analysis_unit": unit_type,
                "namespace": namespace,
                "term_id": term_id,
                "term_name": term_name,
                "background_n": len(frame),
                "foreground_n": foreground_total,
                "control_n": control_total,
                "term_background_n": a + c,
                "term_foreground_n": a,
                "term_control_n": c,
                "foreground_term_fraction": a / foreground_total if foreground_total else np.nan,
                "control_term_fraction": c / control_total if control_total else np.nan,
                "marginal_fraction_difference": (a / foreground_total - c / control_total)
                if foreground_total and control_total
                else np.nan,
                **fisher_exact(a, b, c, d),
                **logistic(frame),
            }
            results.append(row)
    result = pd.DataFrame(results)
    if len(result):
        result["fisher_q"] = np.nan
        result["adjusted_q"] = np.nan
        for (_, namespace), idx in result.groupby(["study", "namespace"], sort=False).groups.items():
            index = np.asarray(list(idx))
            fisher_p = result.loc[index, "fisher_p"].fillna(1.0).astype(float)
            result.loc[index, "fisher_q"] = multipletests(fisher_p, method="fdr_bh")[1]
            adjusted_p = result.loc[index, "adjusted_p"]
            q_all = multipletests(adjusted_p.fillna(1.0).astype(float), method="fdr_bh")[1]
            valid = adjusted_p.notna().to_numpy()
            result.loc[index[valid], "adjusted_q"] = q_all[valid]
        result.sort_values(["namespace", "term_id"], inplace=True)
    return result, pd.DataFrame(terms_all), pd.DataFrame(namespace_rows)

def normalize_concept(term_id: object) -> str:
    value = str(term_id)
    if value.startswith("GO:"):
        return value
    digits = "".join(c for c in value if c.isdigit())
    return digits[-5:] if len(digits) >= 5 else value

def build_display(
    template_path: Path,
    all_results: pd.DataFrame,
    eligibility: pd.DataFrame,
) -> pd.DataFrame:
    template = pd.read_csv(template_path, sep="\t", dtype=str, keep_default_na=False)
    template = template[["study", "species", "ontology", "namespace", "term_id", "term_name"]].copy()
    template.insert(0, "display_order", np.arange(1, len(template) + 1))
    template["concept_id"] = template["term_id"].map(normalize_concept)
    result = all_results.copy()
    result["concept_id"] = result["term_id"].map(normalize_concept)
    columns = [
        "study",
        "namespace",
        "concept_id",
        "background_n",
        "foreground_n",
        "control_n",
        "term_background_n",
        "term_foreground_n",
        "term_control_n",
        "foreground_term_fraction",
        "control_term_fraction",
        "marginal_fraction_difference",
        "fisher_or",
        "fisher_ci_low",
        "fisher_ci_high",
        "fisher_p",
        "fisher_q",
        "adjusted_or",
        "adjusted_ci_low",
        "adjusted_ci_high",
        "adjusted_p",
        "adjusted_q",
        "adjusted_status",
    ]
    merged = template.merge(result[columns], on=["study", "namespace", "concept_id"], how="left")
    eligible = eligibility.copy()
    eligible["concept_id"] = eligible["term_id"].map(normalize_concept)
    merged = merged.merge(
        eligible[["study", "namespace", "concept_id", "term_background_units", "eligibility_status"]].drop_duplicates(
            ["study", "namespace", "concept_id"]
        ),
        on=["study", "namespace", "concept_id"],
        how="left",
    )
    merged["eligibility_status"] = merged["eligibility_status"].replace("", np.nan).fillna(
        "term_absent_from_namespace_membership"
    )
    merged["significant"] = merged["adjusted_q"].lt(0.05)
    merged["adjusted_log2OR"] = np.log2(pd.to_numeric(merged["adjusted_or"], errors="coerce"))
    merged["adjusted_OR"] = merged["adjusted_or"]
    merged["estimability_status"] = merged["adjusted_status"].fillna(merged["eligibility_status"])
    merged["na_reason"] = np.where(
        merged["adjusted_or"].notna(),
        "",
        np.where(
            merged["eligibility_status"].ne("eligible"),
            merged["eligibility_status"],
            merged["estimability_status"],
        ),
    )
    merged.sort_values("display_order", inplace=True)
    if len(merged) != 208:
        raise ValueError(f"fixed display must contain 208 cells, got {len(merged)}")
    return merged
