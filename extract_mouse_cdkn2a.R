#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Seurat))

mouse <- readRDS("D:/IA_AAA_analysis/results/singlecell/mouse_results.rds")
cat("Names:", names(mouse), "\n")
mc <- mouse$mouse_cells
cat("mouse_cells class:", class(mc), "\n")
if (inherits(mc, "Seurat")) {
  cat("Cells:", ncol(mc), "\n")
  cat("Meta cols:", colnames(mc@meta.data), "\n")
  cat("Conditions:\n"); print(table(mc$condition))
  cat("Cell types:\n"); print(table(mc$celltype))
  
  m_expr <- mc[["RNA"]]$data
  cdkn2a_m <- grep("^Cdkn2a$", rownames(mc), ignore.case=TRUE)
  cat("Cdkn2a rows:", length(cdkn2a_m), "\n")
  
  if (length(cdkn2a_m) > 0) {
    cv <- as.numeric(m_expr[cdkn2a_m[1], ])
    # Only Formed condition
    formed <- grepl("formed", mc$condition, ignore.case=TRUE)
    cat("Formed cells:", sum(formed), "\n")
    
    results <- data.frame()
    for (ct in unique(mc$celltype[formed])) {
      mask <- formed & (mc$celltype == ct)
      n <- sum(mask)
      if (n > 5) {
        results <- rbind(results, data.frame(
          species="Mouse", celltype=ct, n=n,
          mean_log1p=mean(cv[mask]),
          pct_expressing=mean(cv[mask] > 0)*100))
      }
    }
    print(results)
    
    # Append to human results
    human_res <- read.delim("D:/IA_AAA_analysis/results/singlecell/cdkn2a_by_celltype_new.tsv")
    final <- rbind(human_res, results)
    write.table(final, "D:/IA_AAA_analysis/results/singlecell/cdkn2a_by_celltype_new.tsv",
                sep="\t", row.names=FALSE, quote=FALSE)
    cat("Saved final results\n")
  }
}
