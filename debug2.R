library(coloc)
library(data.table)

OUT_DIR <- "D:/IA_AAA_analysis/results/coloc"

# Load 9p21.3 data
ia <- fread(file.path(OUT_DIR, "region_chr9_21600000_22600000_IA.tsv"))
eqtl <- fread(file.path(OUT_DIR, "eqtl_Artery_Aorta_chr9_21600000_22600000.tsv"))

ia <- ia[!is.na(beta) & !is.na(se) & se > 0 & !is.na(maf) & maf > 0 & maf < 0.5, ]
eqtl <- eqtl[!is.na(slope) & !is.na(slope_se) & slope_se > 0 & !is.na(maf) & maf > 0 & maf < 0.5, ]

ia$snp_id <- as.character(ia$pos)
eqtl$snp_id <- as.character(eqtl$pos_hg19)

# Pick one gene
gene <- "ENSG00000099810.21"
eqtl_gene <- eqtl[phenotype_id == gene, ]

cat("IA rows:", nrow(ia), "\n")
cat("eQTL gene rows:", nrow(eqtl_gene), "\n")

# Build datasets (like debug script that worked)
common <- intersect(ia$snp_id, eqtl_gene$snp_id)
idx1 <- match(common, ia$snp_id)
idx2 <- match(common, eqtl_gene$snp_id)

eqtl_N <- median(eqtl_gene$ma_count / (2 * eqtl_gene$af), na.rm=TRUE)
cat("eqtl_N:", eqtl_N, "class:", class(eqtl_N), "\n")

ds_ia <- list(
  beta = ia$beta[idx1],
  varbeta = ia$se[idx1]^2,
  N = median(ia$n),
  MAF = ia$maf[idx1],
  type = "cc",
  s = 0.33,
  snp = common
)

ds_eqtl <- list(
  beta = eqtl_gene$slope[idx2],
  varbeta = eqtl_gene$slope_se[idx2]^2,
  N = eqtl_N,
  MAF = eqtl_gene$maf[idx2],
  type = "quant",
  snp = common
)

cat("\nds_ia names:", names(ds_ia), "\n")
cat("ds_eqtl names:", names(ds_eqtl), "\n")
cat("ds_eqtl has s?:", !is.null(ds_eqtl$s), "\n")

cat("\nRunning coloc.abf directly...\n")
res <- coloc.abf(ds_ia, ds_eqtl)
print(res$summary)

# Now try with the function approach
cat("\n\n=== Now testing function approach ===\n")

run_coloc_pair <- function(ds1, ds2) {
  common_snps <- intersect(ds1$snp, ds2$snp)
  idx1 <- match(common_snps, ds1$snp)
  idx2 <- match(common_snps, ds2$snp)

  cat("  In function: common_snps:", length(common_snps), "\n")

  ds1_sub <- list(
    beta = ds1$beta[idx1],
    varbeta = ds1$varbeta[idx1],
    N = ds1$N,
    MAF = ds1$MAF[idx1],
    type = ds1$type,
    snp = ds1$snp[idx1]
  )
  if (!is.null(ds1$s)) ds1_sub$s <- ds1$s

  ds2_sub <- list(
    beta = ds2$beta[idx2],
    varbeta = ds2$varbeta[idx2],
    N = ds2$N,
    MAF = ds2$MAF[idx2],
    type = ds2$type,
    snp = ds2$snp[idx2]
  )
  if (!is.null(ds2$s)) ds2_sub$s <- ds2_sub$s  # BUG? should be ds2$s

  cat("  ds1_sub names:", names(ds1_sub), "\n")
  cat("  ds2_sub names:", names(ds2_sub), "\n")
  cat("  ds2_sub type:", ds2_sub$type, "\n")
  cat("  ds2_sub s:", ds2_sub$s, "\n")

  coloc.abf(ds1_sub, ds2_sub)
}

# Build full-size datasets (like in main script)
ds_ia_full <- list(
  beta = ia$beta,
  varbeta = ia$se^2,
  N = median(ia$n),
  MAF = ia$maf,
  type = "cc",
  s = 0.33,
  snp = ia$snp_id
)

ds_eqtl_full <- list(
  beta = eqtl_gene$slope,
  varbeta = eqtl_gene$slope_se^2,
  N = eqtl_N,
  MAF = eqtl_gene$maf,
  type = "quant",
  snp = eqtl_gene$snp_id
)

cat("Calling function with full datasets...\n")
res2 <- run_coloc_pair(ds_ia_full, ds_eqtl_full)
print(res2$summary)
