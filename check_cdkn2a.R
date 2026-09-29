# Quick check: is CDKN2A in the pseudobulk matrix?
suppressPackageStartupMessages(library(Seurat))
setwd("D:/IA_AAA_analysis")

# Read the pseudobulk DE table
res <- read.table("results/singlecell/pseudobulk_de_synthetic_vs_contractile.tsv",
                  sep="\t", header=TRUE, stringsAsFactors=FALSE)
cat("Total genes in DE table:", nrow(res), "\n")
cat("CDKN2A in table:", "CDKN2A" %in% res$gene, "\n")

# Check if CDKN2A is near the bottom
cdk_rows <- grep("CDKN", res$gene, value=TRUE)
cat("CDKN genes in table:", paste(cdk_rows, collapse=", "), "\n")

# Check the RDS for more info
hr <- readRDS("results/singlecell/human_results.rds")
str(hr)
