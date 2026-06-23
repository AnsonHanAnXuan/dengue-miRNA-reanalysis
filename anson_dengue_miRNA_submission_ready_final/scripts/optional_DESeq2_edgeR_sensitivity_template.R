# Optional external robustness check for GSE150623.
# This script is provided for a supervisor or collaborator with R/Bioconductor installed.
# It was not executed inside the current container because R/Bioconductor are not installed here.
# The Python package already contains a count-aware sensitivity analysis using offset models.

# install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "edgeR"))

library(DESeq2)
library(edgeR)

counts <- read.csv("data_processed/GSE150623_raw_counts_collapsed.csv", check.names = FALSE)
rownames(counts) <- counts[[1]]
counts <- counts[, -1]
meta <- read.csv("data_processed/GSE150623_metadata.csv", check.names = FALSE)
counts <- counts[, meta$sample]
meta$severity <- factor(ifelse(meta$severity_binary == "severe", "severe", "nonsevere"), levels = c("nonsevere", "severe"))

# DESeq2
# Filter lowly detected miRNAs; tune threshold transparently in final submission.
keep <- rowSums(counts >= 5) >= 3
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts[keep, ])), colData = meta, design = ~ severity)
dds <- DESeq(dds)
res <- results(dds, contrast = c("severity", "severe", "nonsevere"))
res_df <- as.data.frame(res)
res_df$miRNA <- rownames(res_df)
write.csv(res_df, "results/GSE150623_DESeq2_optional_results.csv", row.names = FALSE)

# edgeR quasi-likelihood GLM
y <- DGEList(counts = counts[keep, ], group = meta$severity)
y <- calcNormFactors(y)
design <- model.matrix(~ severity, data = meta)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
edger_df <- topTags(qlf, n = Inf)$table
edger_df$miRNA <- rownames(edger_df)
write.csv(edger_df, "results/GSE150623_edgeR_optional_results.csv", row.names = FALSE)
