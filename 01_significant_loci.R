# ============================================================
# 01_significant_loci.R
# Identify genome-wide significant loci (P < 5e-8) via 500kb distance clumping
# for IA (Bakker 2020 Stage1) and AAA (AAAGen 2023)
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

options(datatable.fread.nThread = 4)

RESULTS_DIR <- "D:/IA_AAA_analysis/results"
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

P_THRESH <- 5e-8
CLUMP_BP <- 500000   # 500 kb clumping window
LOCUS_HALF <- 500000 # locus = lead +/- 500 kb

# ------------------------------------------------------------
# Clumping function: distance-based, greedy by P
# Input: data.table with columns chr, pos, snp, beta, se, p
# Output: lead loci table
# ------------------------------------------------------------
clump_loci <- function(dt, p_thresh = P_THRESH, clump_bp = CLUMP_BP) {
  # only significant
  sig <- dt[p < p_thresh & !is.na(p) & is.finite(p)]
  if (nrow(sig) == 0) {
    warning("No significant SNPs found!")
    return(data.table())
  }
  setorder(sig, p)   # most significant first

  leads <- list()
  remaining <- copy(sig)
  i <- 1
  while (nrow(remaining) > 0) {
    lead <- remaining[1]
    # find all SNPs on same chr within clump_bp of lead pos
    same_chr <- remaining[chr == lead$chr]
    in_window <- same_chr[abs(pos - lead$pos) <= clump_bp]
    lead$n_snps_significant <- nrow(in_window)
    leads[[i]] <- lead
    # remove all in-window SNPs (and the lead itself) from remaining
    remaining <- remaining[!(chr == lead$chr & abs(pos - lead$pos) <= clump_bp)]
    i <- i + 1
  }
  leads_dt <- rbindlist(leads)
  setorder(leads_dt, chr, pos)
  return(leads_dt)
}

# ------------------------------------------------------------
# Define locus boundaries: lead +/- 500kb, but shrink to midpoint
# between adjacent leads on same chr
# ------------------------------------------------------------
define_boundaries <- function(leads_dt, half = LOCUS_HALF) {
  setorder(leads_dt, chr, pos)
  leads_dt[, locus_idx := .I]
  leads_dt[, `:=`(start = pos - half, end = pos + half)]

  # adjust boundaries to midpoint between adjacent loci on same chr
  for (c in unique(leads_dt$chr)) {
    idx <- which(leads_dt$chr == c)
    if (length(idx) == 1) next
    for (j in seq_along(idx)) {
      k <- idx[j]
      if (j > 1) {
        prev_pos <- leads_dt$pos[idx[j-1]]
        mid <- (prev_pos + leads_dt$pos[k]) / 2
        if (mid > leads_dt$start[k]) leads_dt$start[k] <- floor(mid)
      }
      if (j < length(idx)) {
        next_pos <- leads_dt$pos[idx[j+1]]
        mid <- (next_pos + leads_dt$pos[k]) / 2
        if (mid < leads_dt$end[k]) leads_dt$end[k] <- ceiling(mid)
      }
    }
  }
  leads_dt[start < 1, start := 1]
  return(leads_dt)
}

# ============================================================
# IA
# ============================================================
cat("=== Reading IA GWAS ===\n")
ia_file <- "D:/IA_AAA_analysis/data/gwas/ia_2020/IA.GWAS.BakkerMK.2020.sumstats.Stage_1.txt.gz"
ia <- fread(ia_file, sep = "\t", header = TRUE,
            select = c("CHR","BP","SNP","A_EFF","A_NONEFF","BETA","SE","P"))
setnames(ia, c("CHR","BP","SNP","A_EFF","A_NONEFF","BETA","SE","P"),
              c("chr","pos","snp","a_eff","a_noneff","beta","se","p"))
cat("IA total SNPs:", nrow(ia), "\n")
cat("IA significant (P<5e-8):", nrow(ia[p < P_THRESH & !is.na(p)]), "\n")

ia_leads <- clump_loci(ia)
ia_leads <- define_boundaries(ia_leads)
ia_leads[, disease := "IA"]
cat("IA loci identified:", nrow(ia_leads), "\n")

# Save full significant SNPs per locus for downstream
ia_sig <- ia[p < P_THRESH & !is.na(p)]
fwrite(ia_sig, file.path(RESULTS_DIR, "ia_significant_snps.tsv"), sep = "\t")

# Output locus table
ia_out <- ia_leads[, .(disease, chr, start, end, lead_SNP = snp,
                       lead_pos = pos, lead_P = p,
                       lead_BETA = beta, lead_SE = se,
                       n_snps_significant, nearest_gene = NA_character_)]
fwrite(ia_out, file.path(RESULTS_DIR, "ia_significant_loci.tsv"), sep = "\t")
cat("Top 10 IA loci:\n")
print(head(ia_out[order(lead_P), .(chr, lead_pos, lead_SNP, lead_P, n_snps_significant)], 10))

# ============================================================
# AAA (large file; read only needed columns)
# ============================================================
cat("\n=== Reading AAA GWAS (large, column-select) ===\n")
aaa_file <- "D:/IA_AAA_analysis/data/gwas/aaa_2023/meta-AAAgen-final-sumstat.txt.gz"
aaa <- fread(aaa_file, sep = "\t", header = TRUE,
             select = c("chr","pos","ref","alt","Effectsize","Effectsize_SD","pvalue"))
setnames(aaa, c("Effectsize","Effectsize_SD","pvalue"), c("beta","se","p"))
# pseudo-SNP id
aaa[, snp := paste0("chr", chr, ":", pos, ":", ref, ":", alt)]
cat("AAA total SNPs:", nrow(aaa), "\n")
cat("AAA significant (P<5e-8):", nrow(aaa[p < P_THRESH & !is.na(p)]), "\n")

aaa_leads <- clump_loci(aaa)
aaa_leads <- define_boundaries(aaa_leads)
aaa_leads[, disease := "AAA"]
cat("AAA loci identified:", nrow(aaa_leads), "\n")

aaa_sig <- aaa[p < P_THRESH & !is.na(p)]
fwrite(aaa_sig, file.path(RESULTS_DIR, "aaa_significant_snps.tsv"), sep = "\t")

aaa_out <- aaa_leads[, .(disease, chr, start, end, lead_SNP = snp,
                         lead_pos = pos, lead_P = p,
                         lead_BETA = beta, lead_SE = se,
                         n_snps_significant, nearest_gene = NA_character_)]
fwrite(aaa_out, file.path(RESULTS_DIR, "aaa_significant_loci.tsv"), sep = "\t")
cat("Top 10 AAA loci:\n")
print(head(aaa_out[order(lead_P), .(chr, lead_pos, lead_SNP, lead_P, n_snps_significant)], 10))

# ============================================================
# Summary
# ============================================================
cat("\n=== SUMMARY ===\n")
cat("IA loci:", nrow(ia_out), "(expected ~17)\n")
cat("AAA loci:", nrow(aaa_out), "(expected ~121)\n")
cat("\nIA chromosome distribution:\n")
print(table(ia_out$chr))
cat("\nAAA chromosome distribution (top 10):\n")
print(sort(table(aaa_out$chr), decreasing = TRUE)[1:10])
cat("\nDone.\n")
