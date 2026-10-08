from __future__ import annotations

import json, math

from pathlib import Path


import matplotlib as mpl

import matplotlib.pyplot as plt

from matplotlib.patches import FancyBboxPatch, Polygon
from plotting_common import export

COL={'ink':'#202c3b','muted':'#536476','teal':'#18766e','purple':'#7129df','amber':'#aa6804','blue':'#2855dc','line':'#d8e1ea','pale':'#f5f8fb'}

def box(ax,x,y,w,h,face='white',edge=COL['line'],radius=.65,lw=.6):
    p=FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={radius}',facecolor=face,edgecolor=edge,linewidth=lw);ax.add_patch(p);return p

def text(ax,x,y,value,size=8,color=None,weight='normal',ha='left',va='center',**kw):
    return ax.text(x,y,value,fontsize=size,color=color or COL['ink'],fontweight=weight,ha=ha,va=va,**kw)

def field(ax,x,y,w,label,value):
    text(ax,x,y-2.3,label,7.2,color=COL['muted']);box(ax,x,y,w,7.0);text(ax,x+2,y+3.5,value,8)
    ax.add_patch(Polygon([[x+w-4,y+2.8],[x+w-2,y+2.8],[x+w-3,y+4.0]],closed=True,color=COL['muted'],lw=0))

def lanes(rows,gap=2):
    result=[];ends=[]
    for item in sorted(rows,key=lambda r:(int(r.get('start',r.get('position'))),int(r.get('end',r.get('position'))))):
        a=int(item.get('start',item.get('position')));b=int(item.get('end',item.get('position')))
        lane=next((i for i,end in enumerate(ends) if a>end+gap),len(ends))
        if lane==len(ends):ends.append(b)
        else:ends[lane]=b
        result.append((item,a,b,lane))
    return result,len(ends)

def draw(data: Path, output: Path) -> None:
    tracks=json.loads((data/'tracks.json').read_text());meta={'counts':json.loads((data/'counts.json').read_text())}
    assert len(tracks['sites'])==21 and len(tracks['ptms'])==48 and len(tracks['structures'])==21
    assert tracks['references'][0]['length']==335
    mpl.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.6,'savefig.facecolor':'white'})
    fw,fh=177.8,221.0;w=169.8
    fig=plt.figure(figsize=(fw/25.4,fh/25.4))
    fig.set_size_inches(504/72,626.4/72,forward=False)
    def axes(name,top,height):
        ax=fig.add_axes([4/fw,(fh-top-height)/fh,w/fw,height/fh],label=name);ax.set_xlim(0,w);ax.set_ylim(height,0);ax.axis('off');ax.set_gid(name);return ax
    a=axes('a',3,82);b=axes('b',90,128)
    text(a,0,1,'a',12,weight='bold',va='top');text(a,7,1,'Atlas search and resource access',10.5,weight='bold',va='top')
    box(a,0,8,w,73,face=COL['pale'],edge=COL['line'],radius=1.3)
    text(a,3,14,'KlaAtlas',13,weight='bold',color=COL['teal'])
    for x,label in [(52,'Explore'),(76,'Studies'),(102,'Tools'),(126,'Figures'),(150,'About')]:text(a,x,14,label,7.5,color=COL['muted'])
    text(a,3,23,'Explore lysine lactylation at the protein coordinate.',10.1,weight='bold')
    box(a,3,29,w-6,8.5,edge='#94a7b8');text(a,6,33.25,'Gene, protein name or accession',8.1,color=COL['muted'])
    box(a,w-28,29,25,8.5,face=COL['teal'],edge=COL['teal']);text(a,w-15.5,33.25,'Search',8.5,color='white',weight='bold',ha='center')
    field(a,3,44,108,'Organism','All organisms')
    box(a,116,44,w-119,7);text(a,119,47.5,'Batch candidates  →',7.7,weight='bold')
    text(a,3,55,'Example record: GAPDH · Homo sapiens · K194',7.5,color=COL['teal'])
    box(a,3,59,53,6.5,face=COL['teal'],edge=COL['teal']);text(a,5,62.25,'Start experimental shortlist  →',7.5,color='white',weight='bold')
    text(a,64,62.25,'Scope & methods  →',7.5,color=COL['teal']);text(a,w-3,62.25,'www.kla-atlas.com',7.4,color=COL['muted'],ha='right')
    counts=meta['counts'];metrics=[('sites','sites'),('proteins','proteins'),('taxa','taxa'),('studies','publications'),('source_tables','source tables')]
    for i,(key,label) in enumerate(metrics):
        x=3+(w-6)*(i+.5)/5;text(a,x,70,f"{counts[key]:,} {label}",7.8,weight='bold',ha='center')
    a.plot([3,w-3],[73.5,73.5],color=COL['line'],lw=.6)
    text(a,3,77,'Current database',7.4,weight='bold');text(a,w-3,77,f"{counts['sites']:,} sites    {counts['evidence']:,} evidence rows",7.4,ha='right')
    text(b,0,1,'b',12,weight='bold',va='top');text(b,7,1,'Reference-bound sequence and evidence tracks',10.5,weight='bold',va='top')
    text(b,1,8,'Human GAPDH  ·  P04406  ·  335 aa',8.4,weight='bold')
    text(b,1,15.5,'From',7.5,color=COL['muted']);box(b,12,12.9,11,5.3);text(b,17.5,15.55,'1',7.5,ha='center')
    text(b,27,15.5,'To',7.5,color=COL['muted']);box(b,34,12.9,14,5.3);text(b,41,15.55,'335',7.5,ha='center')
    for x,label in [(52,'Apply'),(68,'Reset')]:box(b,x,12.9,14,5.3,face=COL['pale']);text(b,x+7,15.55,label,7.0,ha='center')
    text(b,86,15.5,'Zoom',7.2,color=COL['muted']);b.plot([97,109],[15.5,15.5],color='#a9b8c5',lw=1.1);b.scatter([97],[15.5],s=20,c=COL['teal'],zorder=3);text(b,111,15.5,'1×',7.2)
    box(b,124,12.9,w-125,5.3,face=COL['pale']);text(b,126,15.55,'Other PTMs: All types (48)',7.0)
    for x,label in [(1,'Kla sites'),(28,'Functional features'),(69,'Other PTMs'),(99,'Peptide coverage'),(136,'Structure features')]:
        box(b,x,21.5,2.3,2.3,face=COL['teal'],edge=COL['teal'],radius=.2,lw=.3)
        b.plot([x+.4,x+1,x+1.9],[22.6,23.2,22],color='white',lw=.6)
        text(b,x+3.5,22.65,label,7.0)
    x0,x1=31,w-3
    xp=lambda p:x0+(float(p)-1)/334*(x1-x0)
    for tick in [1,49,96,144,192,240,287,335]:
        x=xp(tick);text(b,x,30,str(tick),7.2,color=COL['muted'],ha='center');b.plot([x,x],[45,111],color=COL['line'],lw=.5,ls=(0,(1.5,3)))
    box(b,x0,33,x1-x0,1.5,face='#eaf5f1',edge=COL['teal'],radius=.5,lw=.65);text(b,1,33.75,'Reference',7.6,color=COL['muted'])
    sites=tracks['sites'];assert [s['position'] for s in sites]==sorted(s['position'] for s in sites)
    # Marker size uses log2(evidence_count + 1).
    assert all(s['evidence_count'] + 1 > 0 for s in sites)
    label_sites={5,27,61,86,107,139,162,186,215,251,334}
    for s in sites:
        radius=.28*min(9,3+math.log2(s['evidence_count']+1));b.scatter([xp(s['position'])],[43],s=(2*radius)**2,c=COL['purple'],edgecolors='white',linewidths=.5,zorder=4)
        if s['position'] in label_sites:text(b,xp(s['position']),38.5,'K'+str(s['position']),7.1,ha='center',color=COL['muted'])
    text(b,1,43,'Kla sites\n(21)',7.6,color=COL['muted'],linespacing=1.15)
    ledger={};geometry=[]
    for key,start_y,label,color,mode in [('features',51.5,'Functional\nfeatures',COL['teal'],'interval'),('ptms',71,'Other PTMs\n(48 annotations)',COL['amber'],'point'),('peptides',91,'Peptide\ncoverage',COL['amber'],'interval')]:
        placed,n=lanes(tracks[key]);ledger[key]={'records':len(placed),'lanes':n};center=start_y+(n-1)
        b.plot([x0,x1],[center,center],color=COL['line'],lw=.7,zorder=0);text(b,1,center,label,7.6,color=COL['muted'],linespacing=1.15)
        for r,start,end,lane in placed:
            y=start_y+lane*2;xx=xp(start)
            if mode=='point':b.scatter([xx],[y],s=9.5,c=color,edgecolors='white',linewidths=.5,zorder=3)
            elif start==end:b.plot([xx,xx],[y-.6,y+.6],color=color,lw=1.0,solid_capstyle='round',zorder=3)
            else:box(b,xx,y-.62,xp(end)-xx,1.24,face=color,edge=color,radius=.28,lw=0)
            geometry.append({'track':key,'start':start,'end':end,'lane':lane})
    y=107.5;b.plot([x0,x1],[y,y],color=COL['line'],lw=.7,zorder=0)
    for r in tracks['structures']:b.scatter([xp(r['position'])],[y],s=12,c=COL['blue'],edgecolors='white',linewidths=.5,zorder=3)
    text(b,1,y,'Structure\nfeatures (21)',7.6,color=COL['muted'],linespacing=1.15)
    for x,width,label in [(1,36,'Copy interval FASTA'),(40,41,'Download interval FASTA'),(84,44,'Copy reproducible range'),(131,38,'Export visible interval')]:
        box(b,x,115.3,width,6.4);text(b,x+width/2,118.5,label,7.0,ha='center')
    text(b,1,125,'Marker size reflects evidence-row support. Track positions follow the displayed reference.',7.0,color=COL['muted'])
    fig.canvas.draw()
    base=output/'Figure_5'
    export(fig,base)
    plt.close(fig)
