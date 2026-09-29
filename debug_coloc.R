library(coloc)
library(data.table)

cat("R version:", R.version.string, "\n")
cat("coloc version:", as.character(packageVersion("coloc")), "\n\n")

# Test 1: both quant with simulated data
cat("=== Test 1: quant vs quant (simulated) ===\n")
set.seed(42)
ds1 <- list(beta=rnorm(100), varbeta=runif(100, 0.01, 0.1), N=500, MAF=runif(100, 0.1, 0.5), type="quant", snp=paste0("snp",1:100))
ds2 <- list(beta=rnorm(100), varbeta=runif(100, 0.01, 0.1), N=500, MAF=runif(100, 0.1, 0.5), type="quant", snp=paste0("snp",1:100))
res <- coloc.abf(ds1, ds2)
print(res$summary)

# Test 2: cc vs quant
cat("\n=== Test 2: cc vs quant (simulated) ===\n")
ds1_cc <- list(beta=rnorm(100), varbeta=runif(100, 0.01, 0.1), N=10000, MAF=runif(100, 0.1, 0.5), type="cc", s=0.3, snp=paste0("snp",1:100))
res2 <- coloc.abf(ds1_cc, ds2)
print(res2$summary)

# Test 3: real data - IA vs eQTL for 9p21.3
cat("\n=== Test 3: Real data - IA vs CDKN2B eQTL ===\n")
ia <- fread("D:/IA_AAA_analysis/results/coloc/region_chr9_21600000_22600000_IA.tsv")
eqtl <- fread("D:/IA_AAA_analysis/results/coloc/eqtl_Artery_Aorta_chr9_21600000_22600000.tsv")

gene <- "ENSG00000099810.21"
eqtl_gene <- eqtl[eqtl$phenotype_id == gene, ]

ia$snp_id <- as.character(ia$pos)
eqtl_gene$snp_id <- as.character(eqtl_gene$pos_hg19)

common <- intersect(ia$snp_id, eqtl_gene$snp_id)
cat("Common SNPs:", length(common), "\n")

idx1 <- match(common, ia$snp_id)
idx2 <- match(common, eqtl_gene$snp_id)

eqtl_N <- median(eqtl_gene$ma_count / (2 * eqtl_gene$af))
cat("eQTL N:", eqtl_N, "\n")

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

cat("Running coloc.abf...\n")
res3 <- coloc.abf(ds_ia, ds_eqtl)
print(res3$summary)

cat("\nAll tests passed!\n")
