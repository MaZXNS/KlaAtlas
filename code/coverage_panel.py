"""Aligned taxon coverage columns using the scales in the submitted figure."""
from __future__ import annotations
import matplotlib
matplotlib.use('Agg')
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import numpy as np
import pandas as pd


def ordered_taxa(table:pd.DataFrame)->pd.DataFrame:
    """Keep exact input counts and stable common ordering across all three columns."""
    required=['taxon_id','species','sites','proteins','publications']
    if any(k not in table for k in required):raise ValueError('Missing taxon/count column')
    if table.taxon_id.duplicated().any():raise ValueError('Taxon rows must be unique')
    for key in ['sites','proteins','publications']:
        values=pd.to_numeric(table[key],errors='raise')
        if not np.isfinite(values).all()or(values<1).any()or not(values==np.floor(values)).all():raise ValueError('Expected positive integer taxon counts')
    return table.sort_values(['sites','proteins','taxon_id'],ascending=[False,False,True],kind='stable').reset_index(drop=True)


def draw_taxon_coverage(figure:Figure,table:pd.DataFrame,positions:list[list[float]]|None=None)->list[Axes]:
    """Render existing counts; per-taxon report identities may overlap across taxa."""
    tax=ordered_taxa(table);positions=positions or [[.22,.39475,.10,.3555],[.36,.39475,.10,.3555],[.50,.39475,.10,.3555]]
    labels=[]
    for row in tax.itertuples():
        parts=str(row.species).split();labels.append((parts[0][0]+'. '+parts[1])if len(parts)>1 else str(row.species))
    axes=[];y=np.arange(len(tax));metrics=[('sites','Sites (log10)','#287fa7'),('proteins','Proteins (log10)','#258975'),('publications','Source reports (log10)','#647985')]
    for i,(column,title,color)in enumerate(metrics):
        axis=figure.add_axes(positions[i]);axes.append(axis);values=tax[column].to_numpy(float)
        for yy in y:axis.axhline(yy,color='#edf1f3',lw=.4,zorder=0)
        if column=='publications':
            # The submitted log10 panel measures each bar from 0.5 and places
            # the 100-report tick at x=290 pt on a 504-pt page.
            baseline = 0.5
            left = figure.get_figwidth() * 72 * positions[i][0]
            right = left + figure.get_figwidth() * 72 * positions[i][2]
            start, end = 253.12, 290.0
            log_per_point = (np.log10(100) - np.log10(baseline)) / (end - start)
            limits = baseline * 10 ** ((np.array([left, right]) - start) * log_per_point)
            axis.set_xscale('log')
            axis.set_xlim(*limits)
            axis.set_xticks([1, 10, 100], ['1', '10', '100'])
            axis.minorticks_off()
            height = 3 * len(y) / (positions[i][3] * figure.get_figheight() * 72)
            axis.barh(y, values - baseline, left=baseline, height=height, color=color, zorder=2)
            for yy, value in zip(y, values):
                axis.annotate(f'{int(value)}', (value, yy), xytext=(4, 0),
                              textcoords='offset points', ha='left', va='center',
                              fontsize=7, color='#34454e', clip_on=False)
        else:
            axis.set_xscale('log');axis.set_xlim(.72,float(values.max())*1.35);axis.hlines(y,1,values,color=color,lw=3.0,zorder=2);axis.scatter(values,y,s=7,color=color,zorder=3,clip_on=False);axis.set_xticks([1,100,10000],['1','100','10k']);axis.minorticks_off()
        axis.set_ylim(len(tax)-.45,-.55);axis.set_yticks(y,labels if i==0 else[]);axis.tick_params(axis='y',length=0,pad=6,labelsize=7);axis.tick_params(axis='x',labelsize=7.0,length=2.5,pad=2);axis.set_title(title,fontsize=7.0,pad=5)
        for spine in ['left','right','top']:axis.spines[spine].set_visible(False)
        axis.spines['bottom'].set_color('#aeb9bf');axis.spines['bottom'].set_linewidth(.55)
    return axes
