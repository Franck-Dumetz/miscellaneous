#!/usr/bin/env Rscript

# ----------------------------
# Convert raw counts to CPM
# ----------------------------

# Load edgeR
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR", ask = FALSE)

library(edgeR)

# ---- Input file ----
input_file <- "your_raw_counts.txt"     # Change this to your file path
output_cpm <- "CPM_matrix.txt"
output_logcpm <- "logCPM_matrix.txt"

# ---- Load raw counts ----
counts <- read.delim(input_file, row.names = 1, check.names = FALSE)
print("Raw counts loaded:")
print(dim(counts))

# ---- Create DGEList ----
dge <- DGEList(counts = counts)

# ---- Normalize (TMM) ----
dge <- calcNormFactors(dge)

# ---- Compute CPM and logCPM ----
cpm_matrix <- cpm(dge, normalized.lib.sizes = TRUE)
logcpm_matrix <- cpm(dge, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

# ---- Write output ----
write.table(cpm_matrix, file = output_cpm, sep = "\t", quote = FALSE)
write.table(logcpm_matrix, file = output_logcpm, sep = "\t", quote = FALSE)

cat("CPM and logCPM matrices written to disk.\n")
