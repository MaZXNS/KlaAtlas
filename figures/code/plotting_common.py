"""Small common exports used by the figure drawing scripts."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def export(fig: plt.Figure, stem: Path) -> None:
    """Write the vector figure."""
    stem.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(stem.with_suffix('.pdf'),facecolor='white')
