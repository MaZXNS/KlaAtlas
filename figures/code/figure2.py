from __future__ import annotations


from pathlib import Path

import matplotlib as mpl

import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize


from matplotlib.ticker import FuncFormatter

import numpy as np

import pandas as pd

from plotting_common import export

BLUE='#2B6F9C'

ORANGE='#D9822B'

GREEN='#2B8C72'

GREY='#6F7880'

PAL=['#F8FBFF','#DCEAF5','#9BC3DF','#4D8FBE','#174F7A']

NA='#ECEFF1'

def heat(ax,cells:pd.DataFrame)->None:
    kg=['2','3','4-5'];gg=['1','2-3','4-7','8-15','16+'];m=np.full((3,5),np.nan)
    for i,k in enumerate(kg):
        for j,g in enumerate(gg):
            z=cells[cells.K_group.astype(str).eq(k)&cells.gap_group.eq(g)]
            if len(z):m[i,j]=100*(1-z.equal_table_threshold_fraction.iloc[0])
    cmap=LinearSegmentedColormap.from_list('v1blue',PAL);cmap.set_bad(NA);norm=Normalize(0,18);ax.imshow(np.ma.masked_invalid(m),cmap=cmap,norm=norm,aspect='auto')
    ax.set_xticks(range(5),gg);ax.set_yticks(range(3),kg);ax.set(xlabel='Nearest other K distance (aa)',ylabel='Peptide K count',title='Localization below 0.90 (%)')
    ax.set_xticks(np.arange(-.5,5,1),minor=True);ax.set_yticks(np.arange(-.5,3,1),minor=True);ax.grid(which='minor',color='white',lw=.65);ax.tick_params(which='minor',bottom=False,left=False)
    for i in range(3):
        for j in range(5):ax.text(j,i,'NA' if not np.isfinite(m[i,j]) else f'{m[i,j]:.1f}',ha='center',va='center',fontsize=7,color='white' if np.isfinite(m[i,j]) and m[i,j]>10 else '#17242B')

def draw(data: Path, output: Path) -> None:
    f3=output;f3.mkdir(parents=True,exist_ok=True)
    cov=pd.read_csv(data/'coverage.tsv.gz',sep='\t');pos=pd.read_csv(data/'positions.tsv',sep='\t');om=pd.read_csv(data/'terminal_tests.tsv',sep='\t');bur=pd.read_csv(data/'clustering.tsv',sep='\t');cur=pd.read_csv(data/'length_profiles.tsv',sep='\t');geo=pd.read_csv(data/'peptide_geometry.tsv',sep='\t');cells=pd.read_csv(data/'localization.tsv',sep='\t',dtype={'K_group':str})
    mpl.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,3,figsize=(177.8/25.4,177/25.4));fig.set_size_inches(504/72,501.12/72,forward=False);fig.subplots_adjust(left=.08,right=.98,bottom=.08,top=.95,wspace=.52,hspace=.52)
    for ax,l in zip(axs.flat,'abcdef'):ax.text(-.20,1.06,l,transform=ax.transAxes,fontsize=11,fontweight='bold')
    ax=axs[0,0];ax.hexbin(cov.n_K_evaluable,cov.reported_Kla_fraction,xscale='log',gridsize=(34,22),norm=LogNorm(),mincnt=1,cmap='Blues',linewidths=0,rasterized=True);ax.xaxis.set_major_formatter(FuncFormatter(lambda v,_:f'{v:g}'));ax.set(xlabel='Eligible K per protein',ylabel='Reported Kla fraction',title='Within-protein catalogue coverage')
    ax=axs[0,1];colors={'actual_shortest':BLUE,'actual_median':GREEN,'theoretical_shortest':ORANGE};labels={'actual_shortest':'Observed shortest','actual_median':'Observed median','theoretical_shortest':'Theoretical comparator'}
    for n,g in cur.groupby('model',sort=False):ax.plot(g.x,g.log2_relative_or,color=colors[n],label=labels[n]);ax.fill_between(g.x,g.ci95_lower,g.ci95_upper,color=colors[n],alpha=.13,lw=0)
    ax.axhline(0,color='.65',lw=.5);ax.set(xlabel='Peptide length (aa)',ylabel='Conditional log2 OR\n(relative to 15 aa)',title='Association with peptide length');ax.legend(loc='lower left',frameon=True,facecolor='white',framealpha=.92,edgecolor='#B8C0C6',handlelength=1.3)
    ax=axs[0,2]
    for L,col in [(15,BLUE),(25,ORANGE)]:
        g=geo[geo.peptide_length.eq(L)];ax.plot(g.relative_N_position,g.log2_relative_or,color=col,label=f'{L} aa');ax.fill_between(g.relative_N_position,g.ci95_lower,g.ci95_upper,color=col,alpha=.14,lw=0)
    ax.axhline(0,color='.65',lw=.5);ax.set(xlabel='Position from N to C',ylabel='Conditional log2 OR\n(vs midpoint)',title='Within-peptide position');ax.legend(loc='upper right',frameon=True,facecolor='white',framealpha=.92,edgecolor='#B8C0C6')
    heat(axs[1,0],cells)
    ax=axs[1,1];z=pos[pos.group.eq('all')].sort_values('bin');x=z.midpoint.to_numpy();ax.axvspan(0,.05,color=ORANGE,alpha=.13,lw=0);ax.axvspan(.95,1,color=ORANGE,alpha=.13,lw=0);ax.fill_between(x,(z.expected_fraction+z.null_low)*100,(z.expected_fraction+z.null_high)*100,color='#D9E1E6',label='95% null');ax.plot(x,z.expected_fraction*100,'--',color=GREY,label='All-K');ax.plot(x,z.observed_fraction*100,color=BLUE,label='Reported');ax.set(xlim=(0,1),xlabel='Relative protein position',ylabel='Mean protein fraction (%)',title=f'Global-discovery position (n={int(z.proteins.iloc[0]):,})')
    allrow=om[om.group.eq('all')].iloc[0];terminal=pd.DataFrame({'label':['N','C','N+C'],'estimate':[z.iloc[0].difference,z.iloc[-1].difference,allrow.terminal_difference],'low':[z.iloc[0].null_low,z.iloc[-1].null_low,allrow.terminal_null_low],'high':[z.iloc[0].null_high,z.iloc[-1].null_high,allrow.terminal_null_high]});data_low=float(min(z.observed_fraction.min(),(z.expected_fraction+z.null_low).min())*100);data_high=float(max(z.observed_fraction.max(),(z.expected_fraction+z.null_high).max())*100);span=data_high-data_low;ax.set_ylim(data_low-.10*span,data_high+.80*span)
    ia=ax.inset_axes([.38,.65,.58,.31],facecolor='white');yy=np.arange(3);vals=terminal.estimate.to_numpy()*100;ia.hlines(yy,terminal.low*100,terminal.high*100,color='#9AA6AE',lw=1.1);ia.scatter(vals,yy,s=10,color=[BLUE,ORANGE,'#7A5AA6'],zorder=3);ia.axvline(0,color='#B8B8B8',lw=.5);ia.set_yticks(yy,[f'{lab} {val:+.2f}' for lab,val in zip(terminal.label,vals)]);ia.invert_yaxis();bound=max(1.6,float(np.nanmax(np.abs(np.r_[terminal.low,terminal.high]))*100*1.15));ia.set_xlim(-bound,bound);ia.set_xticks([-1,0,1] if bound<=1.7 else [-2,0,2]);ia.tick_params(length=2,pad=1,labelsize=7);ia.set_xlabel('Excess (pp)',fontsize=7,labelpad=.5)
    ax=axs[1,2]
    for label,col in [('2-4',BLUE),('5-9',GREEN),('10-19',ORANGE),('20+','#7A5AA6')]:
        g=bur[bur.burden.eq(label)];ax.plot(g.width,g.estimate,'o-',ms=3,color=col,label=label);ax.fill_between(g.width,g.ci95_low,g.ci95_high,color=col,alpha=.13,lw=0)
    top=float(bur.ci95_high.max());ax.axhline(0,color='.65',lw=.5);ax.set(xlim=(8,108),ylim=(-.004,top+.020),xlabel='Window width (aa)',ylabel='Mean excess fraction',title='Clustering across fixed scales');ax.set_xticks([10,20,30,50,100]);ax.legend(loc='upper center',ncol=2,frameon=True,facecolor='white',framealpha=.92,edgecolor='#B8C0C6',handlelength=1.2,columnspacing=.8)
    export(fig,f3/'Figure_2')
