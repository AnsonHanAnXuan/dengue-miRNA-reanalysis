from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parents[1]
RES=BASE/'results'
DATA=BASE/'data_processed'
FIGS=BASE/'figures'
FIGS.mkdir(parents=True,exist_ok=True)
primary=pd.read_csv(RES/'target_miRNA_primary_comparison_summary.csv')
meta_study=pd.read_csv(RES/'severity_meta_analysis_per_study_effects.csv')
meta=pd.read_csv(RES/'severity_meta_analysis_random_effects.csv')
counts=pd.read_csv(DATA/'GSE150623_raw_counts_collapsed.csv')
meta150=pd.read_csv(DATA/'GSE150623_metadata.csv')
count_sens=pd.read_csv(RES/'GSE150623_count_model_sensitivity_targets.csv')
libsum=pd.read_csv(RES/'GSE150623_library_size_summary_by_group.csv')

# Figure 1: primary targets, log2CPM/Welch analysis
mirnas=['hsa-miR-574-5p','hsa-miR-1246']
datasets=['GSE190749','GSE150623']
positions=[]; vals=[]; labels=[]; qs=[]
x=0
for ds in datasets:
    for m in mirnas:
        r=primary[(primary.dataset==ds)&(primary.comparison=='SD_vs_NonSevere')&(primary.miRNA==m)].iloc[0]
        positions.append(x); vals.append(float(r['log2FC_groupA_minus_groupB'])); labels.append(m.replace('hsa-','')); qs.append(float(r['q_bh_welch'])); x+=1
    x+=0.8
fig,ax=plt.subplots(figsize=(7.2,4.3))
ax.axhline(0,linewidth=1)
ax.bar(positions,vals)
ax.set_ylabel('log2FC (severe vs non-severe)')
ax.set_title('Primary targets do not survive FDR-corrected validation')
ax.set_xticks(positions)
ax.set_xticklabels(labels,rotation=25,ha='right')
for pos,val,q in zip(positions,vals,qs):
    y=val+(0.12 if val>=0 else -0.12)
    ax.text(pos,y,f'q={q:.3f}',ha='center',va='bottom' if val>=0 else 'top',fontsize=8)
ax.text(0.5,min(vals)-0.65,'GSE190749',ha='center',va='top',fontsize=9)
ax.text(2.8,min(vals)-0.65,'GSE150623',ha='center',va='top',fontsize=9)
ax.set_ylim(min(vals)-0.95,max(vals)+0.7)
fig.tight_layout(); fig.savefig(FIGS/'figure1_primary_targets_fdr.png',dpi=220); plt.close(fig)

# Figure 2: miR-122-5p forest plot
m='hsa-miR-122-5p'
st=meta_study[(meta_study.miRNA==m)&(meta_study.dataset.isin(['GSE190749','GSE150623']))].copy().sort_values('dataset')
pooled=meta[(meta.miRNA==m)&(meta.model=='Main: GSE190749 + GSE150623')].iloc[0]
rows=[]
for _,r in st.iterrows():
    lo=r['hedges_g']-1.96*r['se_hedges_g']; hi=r['hedges_g']+1.96*r['se_hedges_g']
    rows.append((r['study_label'],float(r['hedges_g']),float(lo),float(hi)))
rows.append(('Random-effects pooled',float(pooled['effect']),float(pooled['ci_low']),float(pooled['ci_high'])))
fig,ax=plt.subplots(figsize=(7.2,3.8))
y=np.arange(len(rows))[::-1]
for yi,(lab,e,lo,hi) in zip(y,rows):
    ax.plot([lo,hi],[yi,yi],linewidth=2)
    ax.plot(e,yi,'o',markersize=6)
    ax.text(-2.2,yi,lab,va='center',fontsize=9)
    ax.text(2.1,yi,f'{e:.2f} [{lo:.2f}, {hi:.2f}]',va='center',ha='right',fontsize=9)
ax.axvline(0,linestyle='--',linewidth=1)
ax.set_xlabel("Hedges' g (higher in severe dengue > 0)")
ax.set_yticks([])
ax.set_xlim(-2.3,2.2); ax.set_ylim(-0.6,len(rows)-0.4)
ax.set_title('miR-122-5p shows the most consistent severity-associated signal')
ax.text(-2.2,-0.42,f"Pooled p={pooled['p']:.3f}; I^2={pooled['I2']:.1f}%",fontsize=9)
fig.tight_layout(); fig.savefig(FIGS/'figure2_mir122_forest.png',dpi=220); plt.close(fig)

# Figure 3: miR-1246 detection issue
mi_col=counts.columns[0]
row=counts[counts[mi_col]=='hsa-miR-1246'].iloc[0]
vals=[]
for _,r in meta150.iterrows():
    vals.append({'sample':r['sample'],'clinical_group':r['clinical_group'],'severity_binary':r['severity_binary'],'raw_count':float(row[r['sample']])})
df=pd.DataFrame(vals)
order=['DI','DWS','DS']
fig,ax=plt.subplots(figsize=(6.4,3.8))
for i,g in enumerate(order):
    yy=df[df.clinical_group==g]['raw_count'].values
    jitter=np.linspace(-0.08,0.08,len(yy)) if len(yy)>1 else np.array([0])
    ax.scatter(np.full(len(yy),i)+jitter,np.log10(yy+1),s=28)
ax.set_xticks(range(len(order))); ax.set_xticklabels(['DI','DWS','Severe dengue'])
ax.set_ylabel('log10(raw count + 1)')
ax.set_title('miR-1246 falls below detection in severe GSE150623 samples')
ax.text(2,0.15,'16/16 severe samples = 0 counts',ha='center',fontsize=9)
fig.tight_layout(); fig.savefig(FIGS/'figure3_mir1246_detection.png',dpi=220); plt.close(fig)

# Figure 4: library-size imbalance in GSE150623
qc=meta150.copy(); qc['log10_library_total']=np.log10(qc['library_total'].clip(lower=1))
order=['DI','DWS','DS']
fig,ax=plt.subplots(figsize=(6.4,3.8))
for i,g in enumerate(order):
    yy=qc[qc.clinical_group==g]['log10_library_total'].values
    jitter=np.linspace(-0.08,0.08,len(yy)) if len(yy)>1 else np.array([0])
    ax.scatter(np.full(len(yy),i)+jitter,yy,s=28)
    if len(yy)>0:
        ax.hlines(np.median(yy),i-0.18,i+0.18,linewidth=2)
ax.set_xticks(range(len(order))); ax.set_xticklabels(['DI','DWS','Severe dengue'])
ax.set_ylabel('log10(total miRNA counts)')
ax.set_title('GSE150623 severe samples have much lower library sizes')
fig.tight_layout(); fig.savefig(FIGS/'figure4_gse150623_library_size_qc.png',dpi=220); plt.close(fig)

# Figure 5: count-model sensitivity for selected GSE150623 targets
selected=['hsa-miR-574-5p','hsa-miR-1246','hsa-miR-122-5p','hsa-miR-320a','hsa-miR-486-5p']
rows=[]
for m in selected:
    r=count_sens[count_sens.miRNA==m]
    if r.empty:
        rows.append((m.replace('hsa-',''),np.nan,np.nan))
    else:
        r=r.iloc[0]
        rows.append((m.replace('hsa-',''),float(r['nb2_offset_log2FC']) if pd.notna(r['nb2_offset_log2FC']) else np.nan,float(r['nb2_offset_p']) if pd.notna(r['nb2_offset_p']) else np.nan))
fig,ax=plt.subplots(figsize=(7.2,3.8))
pos=np.arange(len(rows)); vals=[r[1] for r in rows]
ax.axhline(0,linewidth=1)
ax.bar(pos,[0 if np.isnan(v) else v for v in vals])
ax.set_xticks(pos); ax.set_xticklabels([r[0] for r in rows],rotation=25,ha='right')
ax.set_ylabel('NB2 log2FC with library-size offset')
ax.set_title('GSE150623 count-aware sensitivity: miR-122-5p remains positive', pad=14)
ymax = max([v for v in vals if not np.isnan(v)] + [0]) + 0.8
ymin = min([v for v in vals if not np.isnan(v)] + [0]) - 0.4
ax.set_ylim(ymin, ymax)
for p,(lab,v,pval) in zip(pos,rows):
    if np.isnan(v):
        ax.text(p,0.08,'too sparse',ha='center',va='bottom',fontsize=8,rotation=90)
    else:
        y = v + (0.18 if v>=0 else -0.18)
        ax.text(p,y,f'p={pval:.3g}',ha='center',va='bottom' if v>=0 else 'top',fontsize=8)
fig.tight_layout(); fig.savefig(FIGS/'figure5_count_model_sensitivity.png',dpi=220); plt.close(fig)
print('final figures written to', FIGS)
