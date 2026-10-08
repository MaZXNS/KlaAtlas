"""KlaAtlas analysis routines."""
from __future__ import annotations

import argparse

import hashlib

import json

import sqlite3

import subprocess

import tempfile

from pathlib import Path

import numpy as np

import pandas as pd


CONTROLS = [
    "relative_position", "relative_position_sq", "local_flank_basic_fraction",
    "local_flank_acidic_fraction", "local_flank_hydrophobic_fraction",
    "local_flank_pro_gly_fraction",
]

CATEGORIES = ["Nucleus", "Cytoplasm", "Mitochondrion", "ER", "Golgi", "Cell membrane", "Secreted"]

def zscore(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Standardize continuous model columns on the exact analysis population."""
    result = pd.DataFrame(index=frame.index)
    for column in columns:
        values = frame[column].astype(float)
        sd = float(values.std(ddof=0))
        result[column] = (values - float(values.mean())) / sd if sd > 0 else 0.0
    return result

def conditional_fit(frame: pd.DataFrame, predictors: list[str], label: str) -> pd.DataFrame:
    """Fit the fixed Efron conditional logistic specification through R survival::clogit."""
    variation = frame.groupby("protein_unit_id")["is_kla"].agg(["min", "max"])
    keep = variation.index[variation["min"].ne(variation["max"])]
    data = frame.loc[frame["protein_unit_id"].isin(keep)].dropna(subset=predictors + ["is_kla"]).copy()
    model = pd.concat([data[["is_kla", "protein_unit_id"]].reset_index(drop=True), zscore(data, predictors).reset_index(drop=True)], axis=1)
    with tempfile.TemporaryDirectory(prefix="kla_refresh_clogit_") as name:
        directory = Path(name); inp = directory / "input.tsv"; out = directory / "output.tsv"
        model.to_csv(inp, sep="\t", index=False)
        formula = "is_kla ~ " + " + ".join(predictors) + " + strata(protein_unit_id)"
        code = "\n".join([
            "suppressPackageStartupMessages(library(survival))",
            f'd <- read.delim("{inp}", check.names=FALSE)',
            f'fit <- clogit({formula}, data=d, method="efron", control=coxph.control(iter.max=100))',
            's <- summary(fit)$coefficients', 'ci <- confint(fit)',
            'o <- data.frame(predictor=rownames(s),log_or_per_sd=s[,"coef"],se=s[,"se(coef)"],or_per_sd=exp(s[,"coef"]),ci_low=exp(ci[,1]),ci_high=exp(ci[,2]),p_value=s[,"Pr(>|z|)"],log_likelihood=fit$loglik[2])',
            f'write.table(o,file="{out}",sep="\\t",quote=FALSE,row.names=FALSE)',
        ])
        process = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=False)
        if process.returncode:
            raise RuntimeError(process.stderr)
        result = pd.read_csv(out, sep="\t")
    result.insert(0, "model", label)
    result["n_K"] = len(data); result["n_Kla"] = int(data["is_kla"].sum()); result["n_proteins"] = data["protein_unit_id"].nunique()
    result["conditional_method"] = "survival::clogit_Efron_tie_approximation"
    return result

def protein_balanced_difference(frame: pd.DataFrame, replicates: int, seed: int) -> pd.DataFrame:
    """Return the equal-protein RSA difference and protein bootstrap interval."""
    paired = frame.groupby(["protein_unit_id", "is_kla"])["residue_rsa"].mean().unstack().dropna()
    differences = paired[1] - paired[0]
    rng = np.random.default_rng(seed)
    draws = np.array([rng.choice(differences.to_numpy(), len(differences), replace=True).mean() for _ in range(replicates)])
    return pd.DataFrame([{"feature": "residue_rsa", "protein_balanced_difference": differences.mean(),
                          "ci_low": np.quantile(draws, 0.025), "ci_high": np.quantile(draws, 0.975),
                          "n_paired_proteins": len(differences)}])
