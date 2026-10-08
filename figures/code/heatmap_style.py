"""Fixed enrichment display palette and study order."""
import textwrap
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
BLUE = "#2B6F9C"

DARK = "#27323A"

RED = "#C9544D"

HEAT_CMAP = LinearSegmentedColormap.from_list("kla_diverging", [BLUE, "#F7F7F5", RED], N=256)

STUDIES = [
    ("Human hippocampi", "Human_hippocampi", "Hippocampi", "Human"),
    ("OSCC normoxia/hypoxia", "OSCC_normoxia_hypoxia", "OSCC", "Human"),
    ("Acute lymphoblastic leukemia", "Acute_lymphoblastic_leukemia", "ALL", "Human"),
    ("SLC4A7 lung adenocarcinoma", "SLC4A7_lung_adenocarcinoma", "SLC4A7 LUAD", "Human"),
    ("HNRNPC pancreatic cancer", "HNRNPC_pancreatic_cancer", "HNRNPC PDAC", "Human"),
    ("Herpesvirus infection", "Herpesvirus_infection", "Herpesvirus", "Human"),
    ("HLA-F trophoblast", "HLA-F_trophoblast", "HLA-F", "Human"),
    ("Gastrointestinal cancers", "Gastrointestinal_cancers", "GI cancers", "Human"),
    ("PCOS granulosa cells", "PCOS_granulosa_cells", "PCOS", "Human"),
    ("Hypertrophic scar", "Hypertrophic_scar", "Scar", "Human"),
    ("MM.1S WT–LenR", "MM.1S_WT_LenR", "MM.1S", "Human"),
    ("Kawasaki mouse heart", "Kawasaki_mouse_heart", "Kawasaki", "Mouse"),
    ("Serpina3k cardiac ischemia", "Serpina3k_cardiac_ischemia", "Serpina3k", "Mouse"),
]

def panel_label(ax: plt.Axes, label: str) -> None:
    """Place one lower-case panel label with a fixed visual anchor."""
    ax.text(-0.03, 1.025, label, transform=ax.transAxes, fontsize=10.5, fontweight="bold", ha="right", va="bottom")

def wrap_label(label: str, width: int = 29) -> str:
    """Wrap a long term label without altering its wording."""
    return "\n".join(textwrap.wrap(str(label), width=width, break_long_words=False))

HEAT_CMAP.set_bad("#E4E7E9")
plt.rcParams.update({"axes.spines.top":False,"axes.spines.right":False})
