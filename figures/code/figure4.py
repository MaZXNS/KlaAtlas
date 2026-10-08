from __future__ import annotations






from pathlib import Path

import matplotlib as mpl

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd

from matplotlib.colors import Normalize

from matplotlib.patches import Rectangle

from plotting_common import export

def draw_heatmap(ax, cells: pd.DataFrame, selected: pd.DataFrame, recurrent, panel: str, show_namespace: bool) -> None:
    """Reuse the prior heatmap semantics, drawing NA crosses as vector markers."""
    recurrent.panel_label(ax, panel)
    studies = [study for study, _, _, _ in recurrent.STUDIES]
    abbreviations = [abbreviation for _, _, abbreviation, _ in recurrent.STUDIES]
    term_ids = selected["term_id"].astype(str).tolist()
    matrix = np.full((len(term_ids), len(studies)), np.nan)
    significance = np.zeros_like(matrix, dtype=bool)
    for i, term_id in enumerate(term_ids):
        for j, study in enumerate(studies):
            row = cells.loc[(cells["term_id"].astype(str) == term_id) & cells["study"].eq(study)].iloc[0]
            value = pd.to_numeric(pd.Series([row["adjusted_log2OR"]]), errors="coerce").iloc[0]
            if np.isfinite(value):
                matrix[i, j] = np.clip(value, -4, 4)
                significance[i, j] = str(row["significant"]).lower() in {"true", "1"}
    ax.pcolormesh(
        np.arange(len(studies) + 1) - 0.5,
        np.arange(len(term_ids) + 1) - 0.5,
        np.ma.masked_invalid(matrix),
        cmap=recurrent.HEAT_CMAP,
        norm=Normalize(-4, 4),
        shading="flat",
        edgecolors="white",
        linewidth=0.45,
    )
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if not np.isfinite(matrix[i, j]):
                ax.scatter(j, i, marker="x", s=12, c="#8A9196", linewidths=0.75, zorder=4)
            elif significance[i, j]:
                ax.scatter(j, i, s=5.5, c="#111111", linewidths=0, zorder=4)
    labels = [recurrent.wrap_label(name) for name in selected["term_name"]]
    if show_namespace:
        labels = [
            f"[{str(namespace).replace('GO_', '')}] {label}"
            for namespace, label in zip(selected["namespace"], labels)
        ]
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xticks(range(len(studies)), abbreviations, rotation=90, ha="left", rotation_mode="anchor")
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=8)
    ax.tick_params(axis="y", length=0, pad=3)
    for j, (_, _, _, species) in enumerate(recurrent.STUDIES):
        color = "#9ACCC3" if species == "Human" else "#D9B06F"
        ax.add_patch(Rectangle((j - 0.46, -0.78), 0.92, 0.20, facecolor=color, edgecolor="none", clip_on=False))
    for i, row in selected.reset_index(drop=True).iterrows():
        subset = cells.loc[cells["term_id"].astype(str).eq(str(row["term_id"]))]
        values = pd.to_numeric(subset["adjusted_log2OR"], errors="coerce")
        tested = int(subset.loc[values.notna(), "study"].nunique())
        significant_count = int(
            subset.loc[subset["significant"].astype(str).str.lower().isin({"true", "1"}), "study"].nunique()
        )
        ax.text(len(studies) + 0.12, i, f"{significant_count}/{tested}", va="center", ha="left", fontsize=5.6, color=recurrent.DARK)
    if show_namespace:
        for boundary in [2.5, 6.5]:
            ax.axhline(boundary, color="white", lw=2.2, zorder=5)
    ax.set_xlim(-0.5, len(studies) + 1.25)
    ax.set_ylim(len(term_ids) - 0.5, -1.20)
    for spine in ax.spines.values():
        spine.set_visible(False)

def draw(data: Path, output: Path) -> None:
    """Render the non-destructive Figure 4 d/e update."""
    import functional_panels as base
    import heatmap_style as recurrent
    localization_path = data / "localizations.tsv"
    localization = pd.read_csv(localization_path, sep="\t")
    localization["display_label"] = localization["category"].replace({"Unknown": "Unassigned"})
    p057 = pd.read_csv(data/"compartment_effects.tsv", sep="\t")
    absolute = pd.read_csv(data/"rsa_difference.tsv", sep="\t")
    p072 = pd.read_csv(data/"rsa_effects.tsv", sep="\t")
    p072_source = pd.read_csv(data/"rsa_histogram.tsv", sep="\t")
    cells = pd.read_csv(data/"enrichment.tsv", sep="\t", dtype={"term_id": str})
    go_selected = pd.read_csv(data/"go_terms.tsv", sep="\t", dtype={"term_id": str})
    kegg_selected = pd.read_csv(data/"kegg_terms.tsv", sep="\t", dtype={"term_id": str})

    study_order = [study for study, _, _, _ in recurrent.STUDIES]



    go_cells = cells.loc[cells["ontology"].eq("GO")].copy()
    kegg_cells = cells.loc[cells["ontology"].eq("KEGG")].copy()
    go_display = go_selected.copy()
    go_display["term_name"] = [
        "chromatin binding",
        "chromatin",
        "negative regulation of Pol II transcription",
        "ribonucleoprotein complex",
        "spliceosomal complex",
        "mRNA splicing via spliceosome",
        "endoplasmic reticulum membrane",
        "mitochondrion",
    ]
    kegg_display = kegg_selected.copy()

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.6,
            "lines.linewidth": 0.7,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    width_mm, height_mm = 177.8, 200.0
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4), facecolor="white")
    # Preserve the stored physical size across Matplotlib figure-manager rounding.
    fig.set_size_inches(width_mm / 25.4, height_mm / 25.4, forward=False)
    ax_a = fig.add_axes([0.05, 0.73, 0.90, 0.23])
    ax_b = fig.add_axes([0.11, 0.48, 0.33, 0.18])
    ax_c = fig.add_axes([0.62, 0.48, 0.33, 0.18])
    ax_d = fig.add_axes([0.20, 0.060, 0.24, 0.25])
    ax_e = fig.add_axes([0.71, 0.060, 0.24, 0.25])
    base.plot_localization(ax_a, localization)
    base.plot_p057(ax_b, p057)
    for tick in ax_b.get_xticklabels():
        tick.set_rotation(60)
        tick.set_fontsize(7)
        tick.set_ha("right")
        tick.set_rotation_mode("anchor")
    base.plot_p072(ax_c, absolute, p072, p072_source)
    draw_heatmap(ax_d, go_cells, go_display, recurrent, "d", True)
    draw_heatmap(ax_e, kegg_cells, kegg_display, recurrent, "e", False)
    short_studies = ["Hippo", "OSCC", "ALL", "SLC4A7", "HNRNPC", "Herpes", "HLA-F", "GI", "PCOS", "Scar", "MM.1S", "Kawa", "Serp"]
    for axis in (ax_d, ax_e):
        axis.set_xticklabels(short_studies, rotation=90, ha="left", rotation_mode="anchor", fontsize=7)
        axis.tick_params(axis="y", labelsize=7)
        for annotation in axis.texts:
            annotation.set_fontsize(7)
    ax_d.set_yticklabels(
        [
            "[MF] chromatin binding",
            "[CC] chromatin",
            "[BP] negative regulation of\nPol II transcription",
            "[CC] ribonucleoprotein\ncomplex",
            "[CC] spliceosomal complex",
            "[BP] mRNA splicing via\nspliceosome",
            "[CC] endoplasmic reticulum\nmembrane",
            "[CC] mitochondrion",
        ],
        fontsize=7,
    )
    ax_e.set_yticklabels(
        [
            "Spliceosome",
            "Phagocytosis",
            "Coronavirus disease",
            "ATP-dependent chromatin\nremodeling",
            "mRNA surveillance pathway",
            "Lysosome biogenesis",
            "Alzheimer disease",
            "Transcriptional misregulation\nin cancer",
        ],
        fontsize=7,
    )

    colorbar_ax = fig.add_axes([0.35, 0.016, 0.24, 0.010])
    colorbar = mpl.colorbar.ColorbarBase(
        colorbar_ax,
        cmap=recurrent.HEAT_CMAP,
        norm=Normalize(-4, 4),
        orientation="horizontal",
        ticks=[-4, 0, 4],
    )
    colorbar.set_label("Adjusted association with reported Kla, log2 odds ratio", fontsize=7, labelpad=2)
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.ax.tick_params(labelsize=7, length=2, pad=1)
    fig.add_artist(Rectangle((0.055, 0.015), 0.011, 0.010, transform=fig.transFigure, facecolor="#9ACCC3", edgecolor="none"))
    fig.text(0.069, 0.020, "Human", fontsize=7, va="center", color="#39796E")
    fig.add_artist(Rectangle((0.13, 0.015), 0.011, 0.010, transform=fig.transFigure, facecolor="#D9B06F", edgecolor="none"))
    fig.text(0.144, 0.020, "Mouse", fontsize=7, va="center", color="#9A6B25")
    fig.text(0.63, 0.020, "●  adjusted q < 0.05     ×  NA / not estimable", fontsize=7, va="center", color=recurrent.DARK)
    fig.canvas.draw()

    stem = output / "Figure_4"
    export(fig,stem)
    plt.close(fig)
