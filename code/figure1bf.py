from __future__ import annotations



import json





from pathlib import Path


import matplotlib.pyplot as plt


import numpy as np

import pandas as pd


from plotting_common import export
from coverage_panel import draw_taxon_coverage

BLUE = "#287fa7"





TEAL = "#258975"



def render_panels_bf(resource: Path, output: Path, unused: Path|None=None) -> None:
    """Render panels b-f with the three requested gray notes removed."""

    summary = json.loads((resource / "counts.json").read_text())
    tax = pd.read_csv(resource / "taxon_counts.tsv", sep="\t").sort_values("sites", ascending=False).reset_index(drop=True)
    years = pd.read_csv(resource / "publication_years.tsv", sep="\t")
    accumulation = pd.read_csv(resource / "source_accumulation.tsv", sep="\t")
    pair = pd.read_csv(resource / "position_agreement.tsv", sep="\t")
    valid = pair[pair.status.eq("estimable_descriptive")]
    plt.rcParams.update({"font.family":"Arial","font.size":7,"axes.labelsize":7,"axes.titlesize":8,"xtick.labelsize":7,"ytick.labelsize":7,"legend.fontsize":7,"pdf.fonttype":42,"svg.fonttype":"none","axes.linewidth":.6,"axes.spines.top":False,"axes.spines.right":False,"legend.frameon":False})
    fig = plt.figure(figsize=(177.8/25.4, 198/25.4), facecolor="white")
    # Preserve the stored physical size across Matplotlib figure-manager rounding.
    fig.set_size_inches(177.8/25.4, 198/25.4, forward=False)
    def panel_label(axis: plt.Axes, letter: str, title: str, x: float=-.12) -> None:
        axis.text(x,1.04,letter,transform=axis.transAxes,fontweight="bold",fontsize=11,va="bottom")
        axis.text(0,1.04,title,transform=axis.transAxes,fontsize=8,va="bottom")
    b=fig.add_axes([.07,.835,.86,.11]);b.axis("off");panel_label(b,"b","Current discovery atlas",x=-.055)
    metrics=[(summary['sites'],'Sites'),(summary['proteins'],'Protein × taxon'),(summary['publications'],'Source reports'),(summary['source_tables'],'Source tables'),(summary['taxa'],'Taxa'),(summary['evidence'],'Evidence rows')]
    for i,(value,name) in enumerate(metrics):
        x=i/6+.01;b.text(x,.53,f'{value:,}',fontsize=12,fontweight='bold',color=BLUE,ha='left');b.text(x,.28,name,fontsize=7,ha='left')
    c=fig.add_axes([.08,.375,.52,.395]);c.axis("off");panel_label(c,"c","Coverage by taxon",x=-.10)
    insets=draw_taxon_coverage(fig,tax)
    d=fig.add_axes([.66,.405,.27,.34]);panel_label(d,"d","Reports and site coverage",x=-.08);xx=years.year.astype(int);d.bar(xx,years.publications,color='#bad5e2',width=.70);d.set_ylabel('Publications per year');d.set_xlabel('Publication year');d.set_xticks(list(range(int(xx.min()),int(xx.max())+1,2)));d.tick_params(axis='x',rotation=45);dr=d.twinx();dr.plot(xx,years.cumulative_coordinates/1000,'o-',color=BLUE,ms=2.5,lw=1);dr.set_ylabel('Cumulative sites (thousands)',color=BLUE);dr.tick_params(axis='y',colors=BLUE);dr.spines['right'].set_visible(True);d.text(.03,.95,'2026 (partial)',transform=d.transAxes,va='top',fontsize=7)
    e=fig.add_axes([.08,.08,.47,.20]);panel_label(e,"e","Coordinates added by integration",x=-.12);human=accumulation[accumulation.taxon_id.eq(9606)]
    for measure,name,color in [('cumulative_new_protein_coordinates','Newly represented proteins',TEAL),('cumulative_existing_protein_coordinates','Previously represented proteins',BLUE)]:
        z=human[human.measure.eq(measure)];e.plot(z.step,z['median']/1000,label=name,color=color,lw=1);e.fill_between(z.step,z.low/1000,z.high/1000,color=color,alpha=.14,lw=0)
    e.set_xlabel('Human source blocks included');e.set_ylabel('Cumulative coordinates (thousands)');e.legend(loc='upper left',fontsize=7)
    f=fig.add_axes([.66,.08,.27,.20]);panel_label(f,"f","Position agreement across sources",x=-.28);mx=max(valid.observed_distance.max(),valid.expected_distance.max())*1.04;f.scatter(valid.expected_distance,valid.observed_distance,s=4,color=BLUE,alpha=.3,edgecolors='none');f.plot([0,mx],[0,mx],'--',color='#73838b',lw=.6);f.set(xlabel='Expected positional distance',ylabel='Observed positional distance',xlim=(0,mx),ylim=(0,mx));f.text(.96,.06,f'{len(valid):,} pairs\n{len(pair)-len(valid)} unavailable',transform=f.transAxes,ha='right',fontsize=7)
    fig.canvas.draw()
    d.set_ylabel('Source reports per year');d.set_xlabel('Report year')
    # The three coverage axes are comparable; the lower cumulative and
    # distance panels share a row but use intentionally different widths.
    export(fig, output / "Figure_1b_f",axes=[*insets,d,e,f],
           panel_ids=['c_sites','c_proteins','c_reports','d','e','f'],
           row_groups=[['c_sites','c_proteins','c_reports'],['e','f']],column_groups=[])
    plt.close(fig)
