"""KlaAtlas analysis routines."""
from __future__ import annotations

import argparse

from pathlib import Path

import logomaker

import matplotlib

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd

from matplotlib.ticker import MaxNLocator

AA = list("ACDEFGHIKLMNPQRSTVWY")

COLORS = {**{a: "#2F6DB0" for a in "KRH"},
          **{a: "#D94B3D" for a in "DE"},
          **{a: "#26963C" for a in "GSTYCNQ"},
          **{a: "#222222" for a in "AVLIMFWP"}}

KEYS = ["protein_unit_id", "reference_sequence_id", "k_position"]

def preference_heights(counts: pd.DataFrame, background: pd.DataFrame) -> pd.DataFrame:
    """Calculate signed PTM-Logo heights with pseudocount 0 and denominator 0.999999999-q.
    The aligned central K and zero-observation terms have height zero.
    """
    if not counts.index.equals(background.index) or list(counts.columns) != AA or list(background.columns) != AA:
        raise ValueError("Foreground/background must have identical positions and 20 ordered AA columns")
    for matrix in (counts, background):
        if not np.isfinite(matrix.to_numpy()).all() or (matrix < 0).any().any() or (matrix.sum(axis=1) <= 0).any():
            raise ValueError("Counts must be finite, nonnegative, with nonempty rows")
        if 0 not in matrix.index or matrix.loc[0].drop("K").ne(0).any() or matrix.loc[0, "K"] <= 0:
            raise ValueError("The conditioned central position must contain only K")
    p = counts.div(counts.sum(axis=1), axis=0)
    q = background.div(background.sum(axis=1), axis=0)
    flank = counts.index != 0
    pf, qf = p.loc[flank].to_numpy(), q.loc[flank].to_numpy()
    if ((pf > 0) & (qf <= 0)).any() or (qf >= .999999999).any():
        raise ValueError("Flank background has unsupported or degenerate probabilities")
    favored = np.zeros_like(pf)
    disfavored = np.zeros_like(pf)
    observed = pf > 0
    favored[observed] = pf[observed] * np.log2(pf[observed] / qf[observed])
    complement = observed & (pf < 1)
    disfavored[complement] = (1-pf[complement]) * np.log2((1-pf[complement]) / (.999999999-qf[complement]))
    # Equal foreground/background frequencies give zero height.
    values = np.where(pf > qf, np.maximum(favored, 0),
                      np.where(pf < qf, -np.maximum(disfavored, 0), 0))
    result = pd.DataFrame(0., index=counts.index, columns=AA)
    result.loc[flank] = values
    result.index.name = "offset"
    return result

def window_counts(windows: pd.Series) -> pd.DataFrame:
    """Count complete 21-aa K-centered windows."""
    if windows.empty or windows.isna().any() or not windows.str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]{21}").all() or not windows.str[10].eq("K").all():
        raise ValueError("Expected complete canonical 21-aa windows with central K")
    array = np.frombuffer("".join(windows).encode("ascii"), dtype="S1").reshape(-1, 21)
    counts = pd.DataFrame({aa: (array == aa.encode()).sum(axis=0) for aa in AA},
                          index=pd.Index(range(-10, 11), name="offset"))
    return counts

def build_ptm_logo(cases: pd.DataFrame, background_windows: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build the logo using complete K windows from foreground proteins.
    Windows are identified by protein_unit_id, reference ID and K position.
    """
    for table in (cases, background_windows):
        if not set(KEYS + ["window"]).issubset(table.columns):
            raise ValueError("Missing exact site identity or sequence window")
        if table.groupby(KEYS, dropna=False).window.nunique().gt(1).any():
            raise ValueError("Conflicting windows for an exact site identity")
    foreground = cases.drop_duplicates(KEYS).copy()
    counts = window_counts(foreground.window)
    proteins = foreground[KEYS[:2]].drop_duplicates()
    background = background_windows.merge(proteins, on=KEYS[:2], how="inner", validate="many_to_one")
    n_before = len(background)
    background = background.drop_duplicates(KEYS)
    n_unique = len(background)
    background = background.loc[background.window.str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]{21}", na=False) & background.window.str[10].eq("K")].copy()
    binding = foreground[KEYS + ["window"]].merge(background[KEYS + ["window"]], on=KEYS,
                                                  how="left", validate="one_to_one", suffixes=("_fg", "_bg"))
    if binding.window_bg.isna().any() or not binding.window_fg.eq(binding.window_bg).all():
        raise ValueError("Foreground does not bind exactly to current background sequences")
    bg_counts = window_counts(background.window)
    heights = preference_heights(counts, bg_counts)
    stats = pd.DataFrame({"offset": counts.index, "N_foreground": len(foreground),
                          "N_background": len(background), "background_proteins": len(proteins),
                          "background_rows_before_deduplication": n_before,
                          "background_duplicate_rows_removed": n_before-n_unique,
                          "background_incomplete_windows_removed": n_unique-len(background),
                          "favored_bits": heights.clip(lower=0).sum(axis=1).to_numpy(),
                          "disfavored_bits": -heights.clip(upper=0).sum(axis=1).to_numpy()})
    return counts, bg_counts, heights, stats
