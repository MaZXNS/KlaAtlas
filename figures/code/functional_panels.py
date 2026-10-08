from __future__ import annotations








import matplotlib as mpl

import matplotlib.pyplot as plt

from matplotlib.patches import Arc, Circle, ConnectionPatch, Ellipse, PathPatch, Rectangle

from matplotlib.path import Path as MplPath

import numpy as np

import pandas as pd


plt.rcParams['font.family'] = 'sans-serif'

plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']

plt.rcParams['svg.fonttype'] = 'none'

plt.rcParams['pdf.fonttype'] = 42

COLORS = {
    "scar": "#C85C3C",
    "mm1s": "#2878A6",
    "adjusted": "#222222",
    "neutral": "#A9B0B5",
    "unknown": "#D7DADD",
    "nucleus": "#7E68B3",
    "cytoplasm": "#4E9A8C",
    "mitochondrion": "#D98647",
    "er": "#5C8FC1",
    "golgi": "#D2A43C",
    "membrane": "#546E7A",
    "secreted": "#AA6F9E",
    "structural": "#238B8D",
}

def panel_label(ax: mpl.axes.Axes, label: str) -> None:
    """Add an 11 pt bold lowercase panel label."""
    offset = mpl.transforms.ScaledTranslation(-12 / 72, 3 / 72, ax.figure.dpi_scale_trans)
    ax.text(0, 1, label, transform=ax.transAxes + offset, fontsize=11, fontweight="bold", ha="left", va="bottom")

def style_axis(ax: mpl.axes.Axes) -> None:
    """Apply compact journal axis styling."""
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(0.6)
    ax.tick_params(width=0.5, length=2.5, labelsize=7)
    ax.grid(False)

def draw_cell(ax: mpl.axes.Axes) -> None:
    """Draw a detailed editable, non-quantitative cell locator schematic."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    cell = MplPath(
        [(0.09, 0.52), (0.06, 0.76), (0.22, 0.94), (0.45, 0.96),
         (0.68, 0.98), (0.91, 0.83), (0.94, 0.61),
         (0.98, 0.37), (0.85, 0.12), (0.61, 0.07),
         (0.39, 0.03), (0.14, 0.14), (0.09, 0.34),
         (0.075, 0.41), (0.075, 0.47), (0.09, 0.52), (0.09, 0.52)],
        [MplPath.MOVETO] + [MplPath.CURVE4] * 15 + [MplPath.CLOSEPOLY],
    )
    inner = MplPath(
        [(0.115, 0.52), (0.09, 0.75), (0.24, 0.91), (0.46, 0.93),
         (0.67, 0.95), (0.87, 0.81), (0.91, 0.60),
         (0.95, 0.38), (0.82, 0.16), (0.60, 0.105),
         (0.40, 0.07), (0.17, 0.17), (0.12, 0.35),
         (0.105, 0.42), (0.105, 0.47), (0.115, 0.52), (0.115, 0.52)],
        [MplPath.MOVETO] + [MplPath.CURVE4] * 15 + [MplPath.CLOSEPOLY],
    )
    ax.add_patch(PathPatch(cell, facecolor="#F7FAFA", edgecolor=COLORS["membrane"], lw=1.15))
    ax.add_patch(PathPatch(inner, facecolor="none", edgecolor="#9FB1B9", lw=0.55))
    ax.add_patch(Circle((0.43, 0.60), 0.205, facecolor="#EEE9F6", edgecolor=COLORS["nucleus"], lw=0.9))
    ax.add_patch(Circle((0.43, 0.60), 0.177, facecolor="none", edgecolor="#A998CE", lw=0.45))
    ax.add_patch(Circle((0.43, 0.60), 0.060, facecolor="#D4C7E8", edgecolor=COLORS["nucleus"], lw=0.55))
    for angle in np.linspace(0, 2 * np.pi, 10, endpoint=False):
        ax.add_patch(Circle((0.43 + 0.205 * np.cos(angle), 0.60 + 0.205 * np.sin(angle)), 0.008,
                            facecolor="white", edgecolor=COLORS["nucleus"], lw=0.4))
    for x, y, angle in [(0.72, 0.68, 18), (0.70, 0.36, -22)]:
        ax.add_patch(Ellipse((x, y), 0.235, 0.105, angle=angle, facecolor="#F8E1CC", edgecolor=COLORS["mitochondrion"], lw=0.85))
        theta = np.deg2rad(angle)
        for delta in (-0.055, 0.0, 0.055):
            xx, yy = x + delta * np.cos(theta), y + delta * np.sin(theta)
            ax.plot([xx - 0.025 * np.sin(theta), xx + 0.025 * np.sin(theta)],
                    [yy + 0.025 * np.cos(theta), yy - 0.025 * np.cos(theta)], color=COLORS["mitochondrion"], lw=0.55)
    for offset in (-0.065, -0.022, 0.022, 0.065):
        xs = np.linspace(0.14, 0.31, 40)
        ys = 0.44 + offset + 0.018 * np.sin((xs - 0.14) * 35)
        ax.plot(xs, ys, color=COLORS["er"], lw=0.75)
        for xx, yy in zip(xs[::8], ys[::8]):
            ax.add_patch(Circle((xx, yy + 0.012), 0.006, facecolor=COLORS["er"], edgecolor="none"))
    for offset, width in zip((0.00, 0.032, 0.064, 0.096), (0.18, 0.16, 0.135, 0.105)):
        ax.add_patch(Arc((0.67, 0.19 + offset), width, 0.08, theta1=205, theta2=340, color=COLORS["golgi"], lw=1.0))
    for x, y, radius, color in [(0.79, 0.20, 0.018, COLORS["golgi"]), (0.84, 0.27, 0.013, COLORS["golgi"]),
                                 (0.83, 0.52, 0.035, "#DCA0C9"), (0.88, 0.58, 0.020, COLORS["secreted"])]:
        ax.add_patch(Circle((x, y), radius, facecolor="white", edgecolor=color, lw=0.7))
    ax.add_patch(Circle((0.23, 0.72), 0.035, facecolor="#F3D8A1", edgecolor="#B88A30", lw=0.65))
    ax.add_patch(Circle((0.25, 0.77), 0.013, facecolor="#E6B75A", edgecolor="none"))
    ax.plot([0.16, 0.36, 0.59, 0.84], [0.23, 0.15, 0.13, 0.29], color="#83A9A0", lw=0.5)
    ax.plot([0.16, 0.28, 0.56, 0.87], [0.83, 0.88, 0.87, 0.70], color="#83A9A0", lw=0.5)
    ax.plot([0.87, 0.94], [0.59, 0.59], color=COLORS["secreted"], lw=0.6)

def plot_localization(ax: mpl.axes.Axes, frame: pd.DataFrame) -> None:
    """Panel a: cell-to-label connectors aligned with counts and bars."""
    order = ["Nucleus", "Cytoplasm", "Mitochondrion", "ER", "Golgi", "Cell membrane", "Secreted", "Unknown"]
    data = frame.set_index("category").loc[order].reset_index()
    explicit = {
        "Nucleus": COLORS["nucleus"], "Cytoplasm": COLORS["cytoplasm"],
        "Mitochondrion": COLORS["mitochondrion"], "ER": COLORS["er"],
        "Golgi": COLORS["golgi"], "Cell membrane": COLORS["membrane"],
        "Secreted": COLORS["secreted"], "Unknown": COLORS["unknown"],
    }
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("Human protein localization annotations", loc="left", fontsize=8, fontweight="bold", pad=3)
    panel_label(ax, "a")
    ax.text(0.01, 0.93, f"{int(frame['denominator_human_atlas_proteins'].iloc[0]):,} proteins; {int(frame['denominator_human_primary_sites'].iloc[0]):,} reported sites",
            transform=ax.transAxes, fontsize=7, color="#4A5560", ha="left", va="center")
    # Cell diagram.
    cell_ax = ax.inset_axes([0.01, 0.13, 0.28, 0.67])
    draw_cell(cell_ax)
    ax.text(0.15, 0.075, "Cell schematic", transform=ax.transAxes,
            fontsize=7, color="#596269", ha="center", va="center")

    # Localization bars.
    bar_ax = ax.inset_axes([0.69, 0.13, 0.30, 0.65])
    y = np.arange(len(data))[::-1]
    colors = [explicit[name] for name in data["category"]]
    bar_ax.barh(y, data["fraction"] * 100, height=0.62, color=colors, edgecolor="none")
    bar_ax.set_xlim(0, 60)
    bar_ax.set_ylim(-0.7, 7.7)
    bar_ax.set_xticks([0, 25, 50])
    bar_ax.set_yticks([])
    bar_ax.set_xlabel("Proteins annotated (%)", fontsize=7)
    style_axis(bar_ax)

    # One aligned label/count column uses the identical top-to-bottom order.
    ax.text(0.35, 0.84, "Localization", transform=ax.transAxes, fontsize=7,
            fontweight="bold", color="#4A5560", ha="left", va="center")
    ax.text(0.65, 0.865, "proteins /\nlinked Kla sites", transform=ax.transAxes,
            fontsize=7, color="#4A5560", ha="right", va="center", linespacing=1.0)
    anchors = {
        "Nucleus": (0.43, 0.60),
        "Cytoplasm": (0.56, 0.84),
        "Mitochondrion": (0.72, 0.68),
        "ER": (0.23, 0.44),
        "Golgi": (0.68, 0.25),
        "Cell membrane": (0.94, 0.41),
        "Secreted": (0.94, 0.59),
    }
    rails = {
        "Nucleus": 0.303, "Cytoplasm": 0.307, "Mitochondrion": 0.311,
        "ER": 0.315, "Golgi": 0.319, "Cell membrane": 0.323, "Secreted": 0.327,
    }
    for idx, row in data.iterrows():
        yi = y[idx]
        y_frac = 0.13 + 0.65 * ((yi + 0.7) / 8.4)
        category = "Unassigned" if row["category"] == "Unknown" else row["category"]
        color = explicit[row["category"]]
        ax.add_patch(Rectangle((0.35, y_frac - 0.016), 0.018, 0.032,
                               transform=ax.transAxes, facecolor=color, edgecolor="none"))
        ax.text(0.375, y_frac, category, transform=ax.transAxes, fontsize=7,
                color="#202124" if row["category"] != "Unknown" else "#596269",
                fontweight="bold" if row["category"] == "Unknown" else "normal",
                ha="left", va="center", zorder=4)
        ax.text(0.65, y_frac,
                f"{int(row['proteins']):,} / {int(row['sites_on_annotated_proteins']):,}",
                transform=ax.transAxes, fontsize=7, color="#202124", ha="right", va="center", zorder=4)
        if row["category"] != "Unknown":
            anchor = anchors[row["category"]]
            anchor_display = cell_ax.transData.transform(anchor)
            anchor_axes = ax.transAxes.inverted().transform(anchor_display)
            rail = rails[row["category"]]
            # Figure-level ConnectionPatch segments remain above the cell inset
            # and make the complete structure-to-label path visible.
            segments = [
                (anchor, cell_ax.transData, (0.292, anchor_axes[1]), ax.transAxes),
                ((0.292, anchor_axes[1]), ax.transAxes, (rail, anchor_axes[1]), ax.transAxes),
                ((rail, anchor_axes[1]), ax.transAxes, (rail, y_frac), ax.transAxes),
                ((rail, y_frac), ax.transAxes, (0.343, y_frac), ax.transAxes),
            ]
            for xy_a, coords_a, xy_b, coords_b in segments:
                connector = ConnectionPatch(
                    xyA=xy_a, coordsA=coords_a, xyB=xy_b, coordsB=coords_b,
                    arrowstyle="-", color=color, lw=0.55, clip_on=False,
                    zorder=50, shrinkA=0, shrinkB=0,
                )
                ax.figure.add_artist(connector)
            cell_ax.plot(anchor[0], anchor[1], marker="o", ms=2.4,
                         mfc=color, mec="white", mew=0.35, color=color,
                         clip_on=False, zorder=60)

def plot_p057(ax: mpl.axes.Axes, effects: pd.DataFrame) -> None:
    """Panel b: final-scope P057 compartment sequence refit."""
    category_order = ["Nucleus", "Cytoplasm", "Mitochondrion", "ER", "Golgi", "Cell membrane", "Secreted"]
    feature_order = ["kr_fraction", "de_fraction", "hydrophobicity", "within_decile_position"]
    labels = ["K + R fraction", "D + E fraction", "Hydrophobicity", "Within-decile position"]
    matrix = effects.pivot(index="category", columns="feature", values="log2_or").reindex(index=category_order, columns=feature_order)
    qvals = effects.pivot(index="category", columns="feature", values="q_BH28").reindex(index=category_order, columns=feature_order)
    vmax = max(0.8, float(np.nanmax(np.abs(matrix.to_numpy()))))
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(4), labels, fontsize=7, rotation=35, ha="right", rotation_mode="anchor")
    ax.set_yticks(np.arange(7), category_order, fontsize=7)
    ax.tick_params(length=0)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iloc[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:+.2f}" + ("*" if qvals.iloc[i, j] < 0.05 else ""), ha="center", va="center", fontsize=7, color="white" if abs(value) > vmax * 0.55 else "black")
    ax.set_title("Sequence associations by compartment\nGlobal discovery; selected-only reports censored\n*BH28 q < 0.05", loc="left", fontsize=8, fontweight="bold")
    colorbar = ax.figure.colorbar(image, ax=ax, fraction=0.035, pad=0.025)
    colorbar.ax.tick_params(labelsize=7, width=0.4, length=2)
    colorbar.set_label("Association with reported Kla,\nlog2 OR per 1 SD", fontsize=7)
    panel_label(ax, "b")

def plot_p072(ax: mpl.axes.Axes, absolute: pd.DataFrame, effects: pd.DataFrame, source: pd.DataFrame) -> None:
    """Panel c: accepted P072 distributions and summary estimands."""
    ax.axis("off")
    population=effects.loc[effects["model"].str.contains("pLDDT>=90",regex=False)].iloc[0]
    ax.set_title(f"Surface exposure of reported Kla\n{int(population['n_proteins']):,} proteins; {int(population['n_Kla']):,} reported sites", loc="left", fontsize=8, fontweight="bold")
    density = ax.inset_axes([0.03, 0.08, 0.57, 0.40])
    for label, value, color in (("Eligible lysines not\nreported as Kla", 0, "#8094A2"), ("Reported Kla", 1, COLORS["structural"])):
        histogram=source.loc[source['reported_Kla'].eq(value)].sort_values('bin_lower')
        edges=np.r_[histogram.bin_lower.to_numpy(),histogram.bin_upper.iloc[-1]]
        density.stairs(histogram.density.to_numpy(),edges,baseline=0,fill=False,linewidth=1.0,color=color,label=label)
    density.set_xlim(0, 1.25)
    density.set_ylim(0, 2.8)
    density.set_yticks([0, 1, 2])
    density.set_xlabel("Lysine RSA", fontsize=7)
    density.set_ylabel("Density", fontsize=7)
    density.set_title("")
    ax.text(0.03, 0.79, "High-confidence\npredicted structures,\npLDDT ≥90",
            transform=ax.transAxes, fontsize=7, ha="left", va="bottom")
    for y0, edge, label_text in [
        (0.68, "#8094A2", "Eligible lysines not\nreported as Kla"),
        (0.56, COLORS["structural"], "Reported Kla"),
    ]:
        ax.add_patch(Rectangle((0.03, y0), 0.042, 0.045, transform=ax.transAxes,
                               facecolor="none", edgecolor=edge, lw=1.0, clip_on=False))
        ax.text(0.085, y0 + 0.022, label_text, transform=ax.transAxes,
                fontsize=7, ha="left", va="center")
    style_axis(density)
    left = ax.inset_axes([0.69, 0.70, 0.27, 0.08])
    row = absolute.iloc[0]
    left.errorbar([row["protein_balanced_difference"]], [0], xerr=[[row["protein_balanced_difference"] - row["ci_low"]], [row["ci_high"] - row["protein_balanced_difference"]]], fmt="o", color=COLORS["structural"], lw=0.8, ms=4, capsize=2)
    left.axvline(0, color="#777777", lw=0.5)
    left.set_xlim(0, max(0.035, float(row["ci_high"]) * 1.25))
    left.set_yticks([])
    left.set_title(f"Mean RSA difference:\nKla − comparator\n{row['protein_balanced_difference']:.3f} ({row['ci_low']:.3f}, {row['ci_high']:.3f})", fontsize=7)
    style_axis(left)
    left.tick_params(labelsize=7)
    right = ax.inset_axes([0.69, 0.06, 0.27, 0.08])
    rsa = effects.loc[effects["model"].str.contains("pLDDT>=90", regex=False) & effects["predictor"].eq("residue_rsa")].iloc[0]
    right.errorbar([rsa["or_per_sd"]], [0],
                   xerr=[[rsa["or_per_sd"] - rsa["ci_low"]], [rsa["ci_high"] - rsa["or_per_sd"]]],
                   fmt="D", color=COLORS["structural"], lw=0.8, ms=3.4, capsize=1.6)
    right.axvline(1, color="#777777", lw=0.5)
    right.set_xlim(1.065, 1.135)
    right.set_xticks([1.08, 1.12])
    right.set_yticks([])
    right.set_xlabel("Adjusted OR per 1 SD higher RSA", fontsize=7)
    right.set_title(f"{rsa['or_per_sd']:.3f} ({rsa['ci_low']:.3f}, {rsa['ci_high']:.3f})", fontsize=7)
    style_axis(right)
    right.tick_params(labelsize=7)
    panel_label(ax, "c")
