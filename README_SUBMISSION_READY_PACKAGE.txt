Submission-ready draft package: dengue circulating-miRNA public-data reanalysis
================================================================================

Manuscript files
----------------
- manuscript/Anson_Dengue_miRNA_SubmissionReady_Draft.docx
- manuscript/Anson_Dengue_miRNA_SubmissionReady_Draft.pdf

Main scientific framing
-----------------------
This is a public-data reproducibility / negative-validation reanalysis. It should not be described as a clinically validated dengue biomarker paper. The defensible conclusion is:
- miR-574-5p and miR-1246 do not robustly validate across the two main public severity datasets.
- miR-1246 is essentially below detection in GSE150623 and should be interpreted as a platform/detection issue, not a simple biological null.
- GSE307678 is Malaysian warning-sign context only and is not part of the primary severe-dengue validation logic.
- miR-122-5p is the most consistent hypothesis-generating signal, but it still requires prospective validation.

Package contents
----------------
- data_processed/: processed expression matrices and metadata used in the manuscript.
- results/: all result tables used by the manuscript, including FDR summaries, target-miRNA results, meta-analysis tables, and GSE150623 count-aware sensitivity outputs.
- figures/: final manuscript figures.
- scripts/: scripts used for count-model sensitivity, figure generation, DOCX manuscript generation, and QC checks.
- render_qc/: rendered page images and PDF export used for visual layout verification.
- FINAL_QA_REPORT.txt: key analysis and manuscript checks; all checks passed in the final run.

Suggested script order
----------------------
Run from any location with the package preserved as a folder:

1. python scripts/count_model_sensitivity.py
2. python scripts/make_final_figures.py
3. python scripts/create_final_manuscript.py
4. python scripts/final_qc_checks.py

The scripts use paths relative to their own folder, so they do not require hard-coded local paths.

Optional sensitivity analysis
-----------------------------
An optional R template is included:
- scripts/optional_DESeq2_edgeR_sensitivity_template.R

If a supervisor has R/Bioconductor available, adding DESeq2 or edgeR results for GSE150623 would further strengthen the manuscript before submission. The current Python count-aware sensitivity analysis already supports the cautious manuscript conclusion.

Final notes before journal submission
-------------------------------------
Before submitting to any journal, confirm:
- official affiliation wording,
- corresponding author email,
- target journal formatting requirements,
- final funding statement,
- final competing-interest statement,
- whether a mentor/supervisor should be added or acknowledged according to actual contribution.

Do not submit the same manuscript to more than one journal at the same time.
