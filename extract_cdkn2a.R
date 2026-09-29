#!/usr/bin/env Rscript
# Extract CDKN2A expression by cell type from regenerated Seurat objects
# For Figure 5 panels B and C

suppressPackageStartupMessages(library(Seurat))

# --- Human ---
cat("Loading human Seurat object...\n")
human <- readRDS("D:/IA_AAA_analysis/results/singlecell/human_aaa_seurat_regenerated.rds")
cat(sprintf("Human: %d cells\n", ncol(human)))
cat("Cell types:\n"); print(table(human$celltype))
cat("VSMC states:\n"); print(table(human$vsmc_state, human$condition))

# Get CDKN2A expression
expr_data <- human[["RNA"]]$data
cdkn2a_row <- which(rownames(human) == "CDKN2A")
cat(sprintf("CDKN2A found: %d row\n", length(cdkn2a_row)))

if (length(cdkn2a_row) > 0) {
  cdkn2a_vals <- as.numeric(expr_data[cdkn2a_row, ])
} else {
  # Try case-insensitive
  cdkn2a_row <- grep("^CDKN2A$", rownames(human), ignore.case=TRUE)
  cat(sprintf("CDKN2A case-insensitive: %d\n", length(cdkn2a_row)))
  cdkn2a_vals <- as.numeric(expr_data[cdkn2a_row[1], ])
}

# Summary by cell type
celltypes <- unique(human$celltype)
results <- data.frame()
for (ct in celltypes) {
  mask <- human$celltype == ct
  # For VSMC, use vsmc_state
  if (ct == "VSMC") {
    states <- unique(human$vsmc_state[mask])
    for (st in states) {
      smask <- mask & (human$vsmc_state == st)
  n <- sum(smask)
  if (n > 0) {
    mean_expr <- mean(cdkn2a_vals[smask])
    pct_expr <- mean(cdkn2a_vals[smask] > 0) * 100
    results <- rbind(results, data.frame(
      species="Human", celltype=st, n=n,
      mean_log1p=mean_expr, pct_expressing=pct_expr))
  }
    }
  } else {
    n <- sum(mask)
    if (n > 0) {
      mean_expr <- mean(cdkn2a_vals[mask])
      pct_expr <- mean(cdkn2a_vals[mask] > 0) * 100
      results <- rbind(results, data.frame(
        species="Human", celltype=ct, n=n,
        mean_log1p=mean_expr, pct_expressing=pct_expr))
    }
  }
}
cat("\nHuman CDKN2A by cell type:\n")
print(results)

# --- Mouse ---
cat("\nLoading mouse results...\n")
mouse <- readRDS("D:/IA_AAA_analysis/results/singlecell/mouse_results.rds")
cat("Mouse object class:", class(mouse), "\n")
cat("Names:", names(mouse), "\n")

# mouse_results is a list, check structure
if (inherits(mouse, "Seurat")) {
  mouse_obj <- mouse
} else if ("mouse_cells" %in% names(mouse)) {
  mouse_obj <- mouse$mouse_cells
} else {
  # Try first element
  mouse_obj <- mouse[[1]]
}

if (inherits(mouse_obj, "Seurat")) {
  cat(sprintf("Mouse: %d cells\n", ncol(mouse_obj)))
  cat("Meta cols:", colnames(mouse_obj@meta.data), "\n")
  if ("condition" %in% colnames(mouse_obj@meta.data)) {
    cat("Conditions:\n"); print(table(mouse_obj$condition))
  }
  if ("celltype" %in% colnames(mouse_obj@meta.data)) {
    cat("Cell types:\n"); print(table(mouse_obj$celltype))
  }
  
  # Get CDKN2A
  m_expr <- mouse_obj[["RNA"]]$data
  cdkn2a_m <- grep("^Cdkn2a$", rownames(mouse_obj), ignore.case=TRUE)
  cat(sprintf("Cdkn2a rows: %d\n", length(cdkn2a_m)))
  
  if (length(cdkn2a_m) > 0) {
    cdkn2a_m_vals <- as.numeric(m_expr[cdkn2a_m[1], ])
    
    # Only Formed condition (IA aneurysm)
    if ("condition" %in% colnames(mouse_obj@meta.data)) {
      formed_mask <- grepl("formed|IA|aneurysm", mouse_obj$condition, ignore.case=TRUE)
      cat(sprintf("Formed/IA cells: %d\n", sum(formed_mask)))
    } else {
      formed_mask <- rep(TRUE, ncol(mouse_obj))
    }
    
    m_celltypes <- unique(mouse_obj$celltype[formed_mask])
    for (ct in m_celltypes) {
      mask <- formed_mask & (mouse_obj$celltype == ct)
      n <- sum(mask)
      if (n > 5) {
        mean_expr <- mean(cdkn2a_m_vals[mask])
        pct_expr <- mean(cdkn2a_m_vals[mask] > 0) * 100
        results <- rbind(results, data.frame(
          species="Mouse", celltype=ct, n=n,
          mean_log1p=mean_expr, pct_expressing=pct_expr))
      }
    }
  }
} else {
  cat("WARNING: mouse object is not Seurat\n")
}

cat("\n=== FINAL RESULTS ===\n")
print(results)
write.table(results, "D:/IA_AAA_analysis/results/singlecell/cdkn2a_by_celltype_new.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
cat("\nSaved cdkn2a_by_celltype_new.tsv\n")
