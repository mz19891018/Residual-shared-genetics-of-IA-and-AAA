cat(R.version.string, "\n")
pkgs <- c("Seurat", "SeuratObject", "DESeq2", "limma", "edgeR",
          "hdf5r", "Matrix", "ggplot2", "dplyr", "patchwork",
          "data.table", "RColorBrewer", "pheatmap")
for (p in pkgs) {
  cat(sprintf("%-15s: %s\n", p, requireNamespace(p, quietly=TRUE)))
}
