obj <- readRDS('D:/IA_AAA_analysis/results/singlecell/human_results.rds')
cat('Top-level names:', names(obj), '\n\n')

for(nm in names(obj)) {
  x <- obj[[nm]]
  cat('===', nm, '=== class:', class(x), '\n')
  if(inherits(x, 'Seurat')) {
    cat('  cells:', ncol(x), 'features:', nrow(x), '\n')
    cat('  idents:\n'); print(table(Idents(x)))
    cat('  meta cols:', paste(colnames(x@meta.data), collapse=','), '\n')
    if('condition' %in% colnames(x@meta.data)) { cat('  condition:\n'); print(table(x@meta.data$condition)) }
    if('umap' %in% names(x@reductions)) cat('  has UMAP\n')
  } else if(is.data.frame(x) || is.matrix(x)) {
    cat('  dim:', nrow(x), 'x', ncol(x), '\n')
    cat('  head:\n'); print(head(x))
  } else if(is.list(x)) {
    cat('  names:', names(x), '\n')
  } else {
    cat('  length:', length(x), '\n')
    print(head(x))
  }
  cat('\n')
}
