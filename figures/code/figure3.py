from __future__ import annotations







from pathlib import Path

import logomaker


import matplotlib.pyplot as plt

from matplotlib.colors import TwoSlopeNorm

from matplotlib.lines import Line2D

from matplotlib.ticker import MaxNLocator

import numpy as np

import pandas as pd

from plotting_common import export

WIDTH_MM = 177.8

HEIGHT_MM = 220.0

AA = list("ACDEFGHIKLMNPQRSTVWY")

BLUE = "#2B6F9C"

ORANGE = "#D9822B"

GREEN = "#2B8C72"

PURPLE = "#7A5AA6"

GREY = "#6F7880"

FEATURES = ["kr_fraction", "de_fraction", "hydrophobicity", "within_decile_position"]

FEATURE_LABELS = ["Basic (K/R)", "Acidic (D/E)", "Hydrophobicity", "Protein position"]

def panel_label(axis: plt.Axes, letter: str) -> None:
    """Place a consistent panel letter."""
    axis.annotate(
        letter, (0, 1), xycoords="axes fraction", xytext=(-23, 6),
        textcoords="offset points", fontsize=11, fontweight="bold",
        ha="left", va="bottom",
    )

def draw_ptm_logo(axis: plt.Axes, heights: pd.DataFrame, colors: dict[str, str]) -> None:
    """Draw the flank preferences and central K marker."""
    logomaker.Logo(
        heights, ax=axis, color_scheme=colors, font_name="Arial", font_weight="bold",
        stack_order="big_on_top", center_values=False, fade_probabilities=False,
        flip_below=False, shade_below=0, fade_below=0, width=.92, vpad=.01,
        vsep=0, alpha=1, show_spines=False,
    )
    upper = max(float(heights.clip(lower=0).sum(axis=1).max()), .01)
    lower = max(float(-heights.clip(upper=0).sum(axis=1).min()), .01)
    axis.set_xlim(-10.6, 10.6)
    axis.set_ylim(-lower * 1.12, upper * 1.12)
    axis.axhline(0, color="#AAB0B5", lw=.5, zorder=0)
    axis.axvline(0, color="#C2C7CB", lw=.55, zorder=0)
    axis.text(
        0, 1.015, "K", transform=axis.get_xaxis_transform(), ha="center", va="bottom",
        fontsize=9, fontweight="bold", color=colors["K"], clip_on=False,
    )
    axis.set_xticks([-10, -5, 0, 5, 10], ["−10", "−5", "0", "+5", "+10"])
    axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
    axis.set(xlabel="Position relative to Kla-site K", ylabel="Preference (bits)")
    axis.set_title("PTM-Logo sequence preferences around Kla sites", pad=27)

def draw_histone_overlap(
    left: plt.Axes,
    right: plt.Axes,
    effects: pd.DataFrame,
    denominators: pd.DataFrame,
) -> None:
    """Draw only the estimable detectability-overlap sensitivity as panel e."""
    scenario = "detectability_overlap"
    subset = effects.loc[effects.scenario.eq(scenario)].copy()
    if len(subset) != 12 or not subset.estimable.all():
        raise ValueError("Panel e requires the complete 12-row estimable overlap sensitivity")
    population = denominators.loc[denominators.scenario.eq(scenario)].set_index("group")
    limits = (
        min(-.8, float(subset.ci_low.min()) - .25),
        max(1.0, float(subset.ci_high.max()) + .25),
    )
    left.axvline(0, color="#B0B5BA", lw=.6)
    handles: list[Line2D] = []
    for group, color, shift, title in [
        ("verified_histone", PURPLE, -.13, "Histones"),
        ("annotated_other_protein", GREEN, .13, "Other annotated"),
    ]:
        values = subset.loc[subset.group.eq(group)].set_index("feature").loc[FEATURES]
        for i, item in enumerate(values.itertuples()):
            left.errorbar(
                item.log2_or, i + shift,
                xerr=[[item.log2_or - item.ci_low], [item.ci_high - item.log2_or]],
                fmt="o", color=color, ms=3.8, elinewidth=1, capsize=2,
            )
        handles.append(Line2D(
            [], [], marker="o", linestyle="none", color=color,
            label=f'{title} (n={int(population.loc[group, "informative_units"]):,})',
            markersize=4,
        ))
    left.set(
        yticks=range(4), yticklabels=FEATURE_LABELS, ylim=(3.5, -.7), xlim=limits,
        xlabel="Association with reported Kla\n(log2 OR / parent SD)",
    )
    left.text(
        0, 1.23, "Theoretical-peptide-overlap sensitivity",
        transform=left.transAxes, fontsize=8, ha="left", va="bottom",
    )
    left.legend(
        handles=handles, loc="lower left", bbox_to_anchor=(-.01, 1.02),
        ncol=1, fontsize=7, handletextpad=.3, borderaxespad=0,
    )

    differences = subset.loc[subset.group.eq("histone_minus_other")].set_index("feature").loc[FEATURES]
    right.axvline(0, color="#B0B5BA", lw=.6)
    for i, item in enumerate(differences.itertuples()):
        significant = item.q_BH12 < .05
        right.errorbar(
            item.log2_or, i,
            xerr=[[item.log2_or - item.ci_low], [item.ci_high - item.log2_or]],
            fmt="*" if significant else "o", ms=5.5 if significant else 3.8,
            mfc=GREY if significant else "white", mec=GREY, color=GREY,
            elinewidth=1, capsize=2,
        )
    right.set(
        yticks=range(4), ylim=(3.5, -.7), xlim=limits,
        xlabel="Histone − other\n(log2 OR)", title="Group difference",
    )
    right.spines["left"].set_visible(False)
    right.tick_params(axis="y", labelleft=False, length=0)

def draw(data: Path, output: Path) -> None:
    """Validate sources, render the single-page layout, and write provenance files."""
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 7, "axes.labelsize": 7,
        "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 7, "axes.linewidth": .65, "lines.linewidth": .8,
        "patch.linewidth": .5, "pdf.fonttype": 42, "svg.fonttype": "none",
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False,
    })

    from types import SimpleNamespace
    ptm=SimpleNamespace(COLORS={**{a: "#2F6DB0" for a in "KRH"},**{a: "#D94B3D" for a in "DE"},**{a: "#26963C" for a in "GSTYCNQ"},**{a: "#222222" for a in "AVLIMFWP"}})
    heights = pd.read_csv(
        data / "logo.tsv", sep="\t", index_col=0
    ).reindex(columns=AA)
    p51 = pd.read_csv(data / "sequence_effects.tsv", sep="\t")
    p53 = pd.read_csv(data / "symmetry.tsv", sep="\t")
    effects = pd.read_csv(data/"histone_effects.tsv", sep="\t")
    denominators = pd.read_csv(data/"histone_denominators.tsv", sep="\t")
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), facecolor="white")
    # Preserve the stored physical size across Matplotlib figure-manager rounding.
    fig.set_size_inches(WIDTH_MM / 25.4, HEIGHT_MM / 25.4, forward=False)
    outer = fig.add_gridspec(
        3, 1, left=.085, right=.955, bottom=.075, top=.93,
        height_ratios=[1, 1.30, 1.12], hspace=.54,
    )
    axa = fig.add_subplot(outer[0])
    middle = outer[1].subgridspec(1, 2, wspace=.34)
    bgrid = middle[0].subgridspec(1, 2, width_ratios=[1, .045], wspace=.08)
    cgrid = middle[1].subgridspec(1, 2, width_ratios=[1, .045], wspace=.08)
    axb, cb1 = fig.add_subplot(bgrid[0]), fig.add_subplot(bgrid[1])
    axc, cb2 = fig.add_subplot(cgrid[0]), fig.add_subplot(cgrid[1])
    bottom = outer[2].subgridspec(1, 2, width_ratios=[.76, 1.24], wspace=.52)
    axd = fig.add_subplot(bottom[0])
    egrid = bottom[1].subgridspec(1, 2, width_ratios=[.68, .32], wspace=.06)
    axe1 = fig.add_subplot(egrid[0])
    axe2 = fig.add_subplot(egrid[1], sharey=axe1)
    for axis, letter in [(axa, "a"), (axb, "b"), (axc, "c"), (axd, "d"), (axe1, "e")]:
        panel_label(axis, letter)

    draw_ptm_logo(axa, heights, ptm.COLORS)
    offsets = list(range(-10, 0)) + list(range(1, 11))
    matrix = p51.pivot(index="residue", columns="offset", values="log2_or").reindex(index=AA, columns=offsets)
    vmax = float(np.nanmax(np.abs(matrix.to_numpy())))
    image_b = axb.imshow(matrix, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
    q = p51.pivot(index="residue", columns="offset", values="q_bh_400").reindex(index=AA, columns=offsets)
    yy, xx = np.where(q < .05)
    axb.scatter(xx, yy, s=4, c="black", marker=".", linewidths=0)
    axb.set_xticks([0, 4, 9, 10, 15, 19], [-10, -6, -1, 1, 6, 10])
    axb.set_yticks(range(20), AA)
    axb.set(title="Global conditional effects", xlabel="Position relative to Kla-site K", ylabel="Residue")
    colorbar_b = fig.colorbar(image_b, cax=cb1)
    colorbar_b.ax.set_title("log2 OR", fontsize=7, pad=2)

    symmetry = p53.pivot(index="residue", columns="distance", values="right_minus_left").reindex(index=AA, columns=range(1, 11))
    vmax = float(np.nanmax(np.abs(symmetry.to_numpy())))
    image_c = axc.imshow(symmetry, aspect="auto", cmap="PuOr_r", norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax))
    q = p53.pivot(index="residue", columns="distance", values="q_BH200").reindex(index=AA, columns=range(1, 11))
    yy, xx = np.where(q < .05)
    axc.scatter(xx, yy, s=4, c="black", marker=".", linewidths=0)
    axc.set_xticks(range(10), range(1, 11))
    axc.set_yticks(range(20), AA)
    axc.set(title="Global right-minus-left effects", xlabel="Distance from centre", ylabel="Residue")
    colorbar_c = fig.colorbar(image_c, cax=cb2)
    colorbar_c.ax.set_title("Δ log2 OR", fontsize=7, pad=2)

    residue_names = {"K": "Lys (K)", "R": "Arg (R)", "D": "Asp (D)", "E": "Glu (E)"}
    for amino, color in zip("KRDE", [BLUE, "#66A5C8", ORANGE, "#E5A45A"]):
        values = p53.loc[p53.residue.eq(amino)].sort_values("distance")
        axd.plot(values.distance, values.right_minus_left, "o-", ms=2.5, color=color, label=residue_names[amino])
        axd.fill_between(values.distance, values.ci_low, values.ci_high, color=color, alpha=.10, lw=0)
    axd.axhline(0, color="#AAB0B5", lw=.55)
    axd.set(
        xlim=(.7, 10.3), xticks=[1, 3, 5, 7, 9], xlabel="Distance from centre",
        ylabel="Right − left log2 OR",
    )
    axd.legend(ncol=2, loc="lower right", bbox_to_anchor=(1, 1.005), borderaxespad=0)

    draw_histone_overlap(axe1, axe2, effects, denominators)

    export(fig, output/"Figure_3")
    plt.close(fig)
