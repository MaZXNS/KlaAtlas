"""Enrichment heatmap style and run-supplied study order."""
import textwrap
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
BLUE = "#2B6F9C"

DARK = "#27323A"

RED = "#C9544D"

HEAT_CMAP = LinearSegmentedColormap.from_list("kla_diverging", [BLUE, "#F7F7F5", RED], N=256)

STUDIES = []

def set_studies(registry) -> None:
    """Set this run's ordered studies from its validated independent registry."""
    required = {"study", "study_id", "display_label", "species", "study_order"}
    if not required.issubset(registry.columns):
        raise ValueError(f"Study registry missing columns: {sorted(required - set(registry.columns))}")
    ordered = registry.copy()
    ordered["study_order"] = pd.to_numeric(ordered["study_order"], errors="raise")
    ordered = ordered.sort_values("study_order", kind="stable")
    if ordered.empty or ordered.study.duplicated().any():
        raise ValueError("Study registry must contain unique studies")
    global STUDIES
    STUDIES = [
        (str(row.study), str(row.study_id), str(row.display_label), str(row.species))
        for row in ordered.itertuples(index=False)
    ]

def panel_label(ax: plt.Axes, label: str) -> None:
    """Place one lower-case panel label with a fixed visual anchor."""
    ax.text(-0.03, 1.025, label, transform=ax.transAxes, fontsize=10.5, fontweight="bold", ha="right", va="bottom")

def wrap_label(label: str, width: int = 29) -> str:
    """Wrap a long term label."""
    return "\n".join(textwrap.wrap(str(label), width=width, break_long_words=False))

HEAT_CMAP.set_bad("#E4E7E9")
plt.rcParams.update({"axes.spines.top":False,"axes.spines.right":False})
