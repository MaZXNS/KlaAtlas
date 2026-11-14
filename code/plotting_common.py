"""Small common exports used by the figure drawing scripts."""
from pathlib import Path
from typing import Any
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def export(fig: plt.Figure, stem: Path, **alignment_options: Any) -> None:
    """Measure final panel geometry and write editable vectors plus a preview."""
    from figure_alignment import require_matplotlib_panel_alignment
    stem.parent.mkdir(parents=True,exist_ok=True)
    require_matplotlib_panel_alignment(fig,json_out=stem.with_suffix('.alignment.json'),
                                      tolerance_pt=1.5,gutter_tolerance_pt=1.5,
                                      strict=True,**alignment_options)
    fig.savefig(stem.with_suffix('.pdf'),facecolor='white')
    fig.savefig(stem.with_suffix('.svg'),facecolor='white')
    fig.savefig(stem.with_suffix('.png'),dpi=600,facecolor='white')
