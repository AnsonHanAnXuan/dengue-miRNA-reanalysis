from pathlib import Path
import warnings, math
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from statsmodels.discrete.discrete_model import NegativeBinomial

BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'results'
OUT.mkdir(parents=True, exist_ok=True)
counts=pd.read_csv(BASE/'data_processed/GSE150623_raw_counts_collapsed.csv')
meta=pd.read_csv(BASE/'data_processed/GSE150623_metadata.csv')
mi_col=counts.columns[0]
samples=meta['sample'].tolist()
missing=[s for s in samples if s not in counts.columns]
if missing:
    raise SystemExit('Missing samples in count matrix: '+', '.join(missing))
severity=(meta['severity_binary'].str.lower()=='severe').astype(float).values
X_glm=pd.DataFrame({'intercept':1.0, 'severe':severity})
X_nb=sm.add_constant(severity)
offset=np.log(meta['library_total'].astype(float).clip(lower=1).values)

def bh(vals):
    vals=np.array(vals,dtype=float)
    ok=np.isfinite(vals)
    q=np.full_like(vals, np.nan, dtype=float)
    if ok.sum()>0:
        q[ok]=multipletests(vals[ok], method='fdr_bh')[1]
    return q

def fit_poisson_robust(y):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            fit=sm.GLM(y, X_glm, family=sm.families.Poisson(), offset=offset).fit(maxiter=100, disp=0, cov_type='HC0')
        beta=float(fit.params['severe']); se=float(fit.bse['severe']); p=float(fit.pvalues['severe'])
        return beta/math.log(2), se/math.log(2), p, 'ok'
    except Exception as e:
        return np.nan, np.nan, np.nan, 'failed:'+str(e)[:100]

def fit_nb2(y):
    if y.sum() <= 0 or (y > 0).sum() < 2:
        return np.nan, np.nan, np.nan, np.nan, 'too_sparse'
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            fit=NegativeBinomial(y, X_nb, offset=offset).fit(disp=0, maxiter=500)
        beta=float(fit.params[1]); se=float(fit.bse[1]); p=float(fit.pvalues[1]); alpha=float(fit.params[-1]) if len(fit.params)>2 else np.nan
        return beta/math.log(2), se/math.log(2), p, alpha, 'ok'
    except Exception as e:
        return np.nan, np.nan, np.nan, np.nan, 'failed:'+str(e)[:100]

# All-miRNA Poisson robust count-offset model for an all-feature BH sensitivity check.
rows=[]
for _,row in counts.iterrows():
    mirna=row[mi_col]
    y=row[samples].astype(float).values
    total=float(y.sum()); detected=int((y>0).sum())
    if total==0 or detected<2:
        rows.append({'miRNA':mirna,'total_count':total,'detected_samples':detected,'poisson_offset_log2FC':np.nan,'poisson_offset_se':np.nan,'poisson_offset_p':np.nan,'status':'too_sparse'})
        continue
    lfc,se,p,status=fit_poisson_robust(y)
    rows.append({'miRNA':mirna,'total_count':total,'detected_samples':detected,'poisson_offset_log2FC':lfc,'poisson_offset_se':se,'poisson_offset_p':p,'status':status})
allp=pd.DataFrame(rows)
allp['poisson_offset_q_bh']=bh(allp['poisson_offset_p'])
allp.to_csv(OUT/'GSE150623_poisson_offset_sensitivity_all_miRNAs.csv',index=False)

# Target-level NB2 sensitivity with estimated dispersion.
targets=['hsa-miR-574-5p','hsa-miR-1246','hsa-miR-122-5p','hsa-miR-320a','hsa-miR-486-5p','hsa-miR-30d-5p','hsa-miR-424-5p','hsa-miR-1303']
trs=[]
for _,row in counts[counts[mi_col].isin(targets)].iterrows():
    mirna=row[mi_col]
    y=row[samples].astype(float).values
    nb_lfc,nb_se,nb_p,nb_alpha,nb_status=fit_nb2(y)
    pois=allp[allp.miRNA==mirna].iloc[0]
    trs.append({'miRNA':mirna,'total_count':float(y.sum()),'detected_samples':int((y>0).sum()),
                'nb2_offset_log2FC':nb_lfc,'nb2_offset_se':nb_se,'nb2_offset_p':nb_p,'nb2_alpha':nb_alpha,'nb2_status':nb_status,
                'poisson_offset_log2FC':pois.poisson_offset_log2FC,'poisson_offset_se':pois.poisson_offset_se,
                'poisson_offset_p':pois.poisson_offset_p,'poisson_offset_q_bh_all_miRNAs':pois.poisson_offset_q_bh,'poisson_status':pois.status})
targ=pd.DataFrame(trs)
targ['nb2_offset_q_bh_targets_only']=bh(targ['nb2_offset_p'])
targ.to_csv(OUT/'GSE150623_count_model_sensitivity_targets.csv', index=False)

# Sample/library-size QC.
qc=meta.copy(); qc['log10_library_total']=np.log10(qc['library_total'].clip(lower=1))
qc.to_csv(OUT/'GSE150623_sample_library_size_qc.csv', index=False)
summary=qc.groupby(['clinical_group','severity_binary']).agg(n=('sample','count'), median_library_total=('library_total','median'), min_library_total=('library_total','min'), max_library_total=('library_total','max')).reset_index()
summary.to_csv(OUT/'GSE150623_library_size_summary_by_group.csv', index=False)
print('target count model sensitivity')
print(targ.to_string(index=False))
print('\nsummary')
print(summary.to_string(index=False))
