if (!requireNamespace("variancePartition", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("variancePartition")
}

library(variancePartition)
library(edgeR)         # optional for CPM calculation
library(limma)         # needed for voom
library(tidyverse)

# read expression data
expr <- read.delim("PATH/transcript_counts_newgff.txt")

# read gene features
features <- read.delim("PATH/gene_features.tsv")

# merge
dat <- merge(expr, features, by = "gene_id")

dat$PTU_ID <- as.factor(dat$PTU_ID)
dat$position_in_PTU <- as.numeric(dat$position_in_PTU)
dat$gene_length <- as.numeric(dat$gene_length)

# set gene_id as rownames
exprMat <- as.matrix(dat$CPM)
rownames(exprMat) <- dat$gene_id
colnames(exprMat) <- "CPM"
y <- dat$CPM

# Define the model
form <- ~ (1|PTU_ID) + gene_length + position_in_PTU #(1|PTU_ID) → random effect, gene_length → fixed effect, position_in_PTU → fixed effect

# Fit the model
varPart <- fitExtractVarPartModel(
  y,
  form,
  dat
)

# calculate average variance explained
avgVar <- colMeans(varPart)

# Show results
print(varPart)
plotVarPart(varPart)
print(avgVar)
