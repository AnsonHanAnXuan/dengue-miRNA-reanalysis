from pathlib import Path
import json, math, re, zipfile
import pandas as pd
from docx import Document

BASE=Path(__file__).resolve().parents[1]
RES=BASE/'results'
DATA=BASE/'data_processed'
DOC=BASE/'manuscript'/'Anson_Dengue_miRNA_SubmissionReady_Draft.docx'
report=[]

def ok(name, condition, detail=''):
    report.append({'check':name,'status':'PASS' if condition else 'FAIL','detail':detail})

# Files present
required=[
RES/'target_miRNA_primary_comparison_summary.csv',
RES/'severity_meta_analysis_random_effects.csv',
RES/'GSE150623_count_model_sensitivity_targets.csv',
RES/'GSE150623_library_size_summary_by_group.csv',
DATA/'GSE150623_raw_counts_collapsed.csv',
DATA/'GSE150623_metadata.csv',
DOC]
for p in required:
    ok(f'file present: {p.name}', p.exists(), str(p))

primary=pd.read_csv(RES/'target_miRNA_primary_comparison_summary.csv')
allm=pd.read_csv(RES/'all_datasets_all_miRNAs_all_comparisons.csv')
meta=pd.read_csv(RES/'severity_meta_analysis_random_effects.csv')
counts=pd.read_csv(DATA/'GSE150623_raw_counts_collapsed.csv')
meta150=pd.read_csv(DATA/'GSE150623_metadata.csv')
cs=pd.read_csv(RES/'GSE150623_count_model_sensitivity_targets.csv')
libsum=pd.read_csv(RES/'GSE150623_library_size_summary_by_group.csv')

# Sample match
samples=meta150['sample'].tolist()
ok('GSE150623 metadata samples all in count matrix', all(s in counts.columns for s in samples), f'{len(samples)} samples checked')
ok('GSE150623 sample count equals 39', len(samples)==39, f'n={len(samples)}')
ok('GSE150623 severe/non-severe count equals 16/23', ((meta150.severity_binary=='severe').sum()==16 and (meta150.severity_binary=='non-severe').sum()==23), f"severe={(meta150.severity_binary=='severe').sum()}, nonsevere={(meta150.severity_binary=='non-severe').sum()}")

# FDR checks
sub190=allm[(allm.dataset=='GSE190749')&(allm.comparison=='SD_vs_NonSevere')]
ok('GSE190749 tested miRNA count is 798', len(sub190)==798, f'n={len(sub190)}')
for m in ['hsa-miR-574-5p','hsa-miR-1246']:
    r=primary[(primary.dataset=='GSE190749')&(primary.comparison=='SD_vs_NonSevere')&(primary.miRNA==m)].iloc[0]
    ok(f'GSE190749 {m} q=0.941', abs(float(r.q_bh_welch)-0.941307)<1e-4, f'q={r.q_bh_welch}')
sub150=allm[(allm.dataset=='GSE150623')&(allm.comparison=='SD_vs_NonSevere')]
ok('GSE150623 tested miRNA count is 766', len(sub150)==766, f'n={len(sub150)}')
ok('GSE150623 has FDR-significant miRNAs overall', sub150.q_bh_welch.lt(0.05).sum()==150, f"q<0.05={sub150.q_bh_welch.lt(0.05).sum()}")
for m in ['hsa-miR-574-5p','hsa-miR-1246','hsa-miR-122-5p']:
    r=primary[(primary.dataset=='GSE150623')&(primary.comparison=='SD_vs_NonSevere')&(primary.miRNA==m)].iloc[0]
    ok(f'GSE150623 {m} normalized q not significant', float(r.q_bh_welch)>0.05, f'q={r.q_bh_welch}')

# miR-1246 detection
mi_col=counts.columns[0]
row=counts[counts[mi_col]=='hsa-miR-1246'].iloc[0]
sev=meta150[meta150.severity_binary=='severe']['sample'].tolist(); non=meta150[meta150.severity_binary=='non-severe']['sample'].tolist()
sev_det=sum(float(row[s])>0 for s in sev); non_det=sum(float(row[s])>0 for s in non)
ok('miR-1246 detected in zero severe GSE150623 samples', sev_det==0, f'{sev_det}/16')
ok('miR-1246 detected in one non-severe GSE150623 sample', non_det==1, f'{non_det}/23')

# miR-122 meta
mir122=meta[(meta.miRNA=='hsa-miR-122-5p')&(meta.model=='Main: GSE190749 + GSE150623')].iloc[0]
ok('miR-122 pooled effect approximately 0.71', abs(float(mir122.effect)-0.709607)<1e-4, f"g={mir122.effect}")
ok('miR-122 pooled p approximately 0.017', abs(float(mir122.p)-0.017290)<1e-4, f"p={mir122.p}")
ok('miR-122 I2 approximately 10.6', abs(float(mir122.I2)-10.641777)<1e-3, f"I2={mir122.I2}")

# Count sensitivity checks
r122=cs[cs.miRNA=='hsa-miR-122-5p'].iloc[0]
r574=cs[cs.miRNA=='hsa-miR-574-5p'].iloc[0]
r1246=cs[cs.miRNA=='hsa-miR-1246'].iloc[0]
ok('GSE150623 count model miR-122 positive and significant', float(r122.nb2_offset_log2FC)>0 and float(r122.nb2_offset_p)<0.001, f"NB2 log2FC={r122.nb2_offset_log2FC}, p={r122.nb2_offset_p}")
ok('GSE150623 count model miR-574 not significant', float(r574.nb2_offset_p)>0.05, f"p={r574.nb2_offset_p}")
ok('GSE150623 count model miR-1246 too sparse', str(r1246.nb2_status)=='too_sparse', f"status={r1246.nb2_status}")

# Library-size imbalance
lib_ds=libsum[libsum.clinical_group=='DS'].iloc[0]
lib_di=libsum[libsum.clinical_group=='DI'].iloc[0]
lib_dws=libsum[libsum.clinical_group=='DWS'].iloc[0]
ok('GSE150623 severe median library smaller than DI and DWS', lib_ds.median_library_total < lib_di.median_library_total and lib_ds.median_library_total < lib_dws.median_library_total, f"DS={lib_ds.median_library_total}, DI={lib_di.median_library_total}, DWS={lib_dws.median_library_total}")

# Manuscript text checks
text='\n'.join(p.text for p in Document(DOC).paragraphs)
for bad in ['STUDENT WRITES','WHAT TO WRITE','GEO2R Results Summary','[ADD','[confirm','TODO','Lorem ipsum']:
    ok(f'manuscript has no placeholder token: {bad}', bad not in text, '')
ok('manuscript states GSE307678 not primary severe validation', 'GSE307678 was not used as primary severe-dengue validation' in text or 'GSE307678 was excluded from the primary analysis' in text, '')
ok('manuscript includes AI assistance disclosure', 'AI assistance disclosure' in text, '')
ok('manuscript includes miR-1246 detection explanation', 'miR-1246 was essentially undetectable in GSE150623' in text, '')
ok('manuscript includes count-aware sensitivity', 'count-aware sensitivity' in text, '')

fail=[r for r in report if r['status']!='PASS']
report_path=BASE/'FINAL_QA_REPORT.txt'
with report_path.open('w') as f:
    f.write('Final QA report for Anson dengue miRNA manuscript package\n')
    f.write('='*64+'\n\n')
    for r in report:
        f.write(f"[{r['status']}] {r['check']} -- {r['detail']}\n")
    f.write('\nSummary: '+('ALL CHECKS PASSED' if not fail else f'{len(fail)} CHECKS FAILED')+'\n')
json_path=BASE/'results/final_qc_report.json'
json_path.write_text(json.dumps(report, indent=2))
print(report_path)
print('failures',len(fail))
if fail:
    for r in fail: print(r)
    raise SystemExit(1)
