# ============================================================
# 02_locus_overlap.R
# Find overlapping IA-AAA loci, classify Type A/B/C
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(IRanges)
})

RESULTS_DIR <- "D:/IA_AAA_analysis/results"

# Load loci
ia_loci <- fread(file.path(RESULTS_DIR, "ia_significant_loci.tsv"))
aaa_loci <- fread(file.path(RESULTS_DIR, "aaa_significant_loci.tsv"))
cat("IA loci:", nrow(ia_loci), " AAA loci:", nrow(aaa_loci), "\n")

# Load munged sumstats for direction comparison (HapMap3, aligned to A1)
cat("Loading munged sumstats for direction comparison...\n")
ia_m <- fread("D:/IA_AAA_analysis/data/qc/IA_2020.sumstats.gz")
aaa_m <- fread("D:/IA_AAA_analysis/data/qc/AAA_2023.sumstats.gz")
setnames(ia_m, c("BETA","P"), c("IA_BETA_m","IA_P_m"))
setnames(aaa_m, c("BETA","P"), c("AAA_BETA_m","AAA_P_m"))
# keep only SNP, aligned beta
ia_m_dir <- ia_m[, .(SNP, IA_BETA_m, IA_A1 = A1, IA_A2 = A2)]
aaa_m_dir <- aaa_m[, .(SNP, AAA_BETA_m, AAA_A1 = A1, AAA_A2 = A2)]
dir_compare <- merge(ia_m_dir, aaa_m_dir, by = "SNP")
cat("Munged SNPs present in both:", nrow(dir_compare), "\n")

# ------------------------------------------------------------
# Build GRanges and find overlaps
# ------------------------------------------------------------
ia_gr <- GRanges(seqnames = paste0("chr", ia_loci$chr),
                 ranges = IRanges(start = ia_loci$start, end = ia_loci$end),
                 lead_snp = ia_loci$lead_SNP,
                 lead_pos = ia_loci$lead_pos,
                 lead_p = ia_loci$lead_P,
                 lead_beta = ia_loci$lead_BETA)

aaa_gr <- GRanges(seqnames = paste0("chr", aaa_loci$chr),
                  ranges = IRanges(start = aaa_loci$start, end = aaa_loci$end),
                  lead_snp = aaa_loci$lead_SNP,
                  lead_pos = aaa_loci$lead_pos,
                  lead_p = aaa_loci$lead_P,
                  lead_beta = aaa_loci$lead_BETA)

# findOverlaps: which AAA loci overlap which IA loci
hits <- findOverlaps(ia_gr, aaa_gr, type = "any", ignore.strand = TRUE)
cat("Overlapping locus pairs:", length(hits), "\n")

ov <- data.table(
  ia_idx = queryHits(hits),
  aaa_idx = subjectHits(hits)
)
ov <- cbind(ov,
  ia_loci[ov$ia_idx, .(chr = chr, ia_start = start, ia_end = end,
                       ia_leadSNP = lead_SNP, ia_leadPos = lead_pos,
                       ia_leadP = lead_P, ia_leadBETA = lead_BETA,
                       ia_leadSE = lead_SE, ia_nSig = n_snps_significant)],
  aaa_loci[ov$aaa_idx, .(aaa_start = start, aaa_end = end,
                         aaa_leadSNP = lead_SNP, aaa_leadPos = lead_pos,
                         aaa_leadP = lead_P, aaa_leadBETA = lead_BETA,
                         aaa_leadSE = lead_SE, aaa_nSig = n_snps_significant)]
)
# overlap interval
ov[, `:=`(ov_start = pmax(ia_start, aaa_start),
          ov_end   = pmin(ia_end, aaa_end))]
ov[, ov_width := ov_end - ov_start + 1]
ov[, lead_distance := abs(ia_leadPos - aaa_leadPos)]

# ------------------------------------------------------------
# Direction concordance
# Try: look up IA lead rsID in both munged files
# ------------------------------------------------------------
cat("Computing direction concordance...\n")
ov[, direction_concordant := NA_character_]  # "concordant"/"discordant"/"unknown"
ov[, direction_note := NA_character_]

for (i in seq_len(nrow(ov))) {
  ia_snp <- ov$ia_leadSNP[i]
  # only works if IA lead is rsID
  if (!grepl("^rs", ia_snp)) {
    ov$direction_note[i] <- "IA lead not rsID"
    next
  }
  match <- dir_compare[SNP == ia_snp]
  if (nrow(match) == 0) {
    ov$direction_note[i] <- "IA lead not in munged intersection"
    next
  }
  b1 <- match$IA_BETA_m[1]
  b2 <- match$AAA_BETA_m[1]
  if (is.na(b1) || is.na(b2) || b1 == 0 || b2 == 0) {
    ov$direction_note[i] <- "missing beta"
    next
  }
  # Note: munged files align to A1, but A1 might differ between studies.
  # We check alleles: if IA_A1 == AAA_A1, same sign = concordant.
  # If IA_A1 == AAA_A2 (strand flip), opposite sign = concordant.
  a1i <- match$IA_A1[1]; a2i <- match$IA_A2[1]
  a1a <- match$AAA_A1[1]; a2a <- match$AAA_A2[1]
  same_allele <- (a1i == a1a && a2i == a2a)
  flipped <- (a1i == a2a && a2i == a1a)
  if (same_allele) {
    conc <- sign(b1) == sign(b2)
    ov$direction_note[i] <- "alleles aligned A1"
  } else if (flipped) {
    conc <- sign(b1) != sign(b2)
    ov$direction_note[i] <- "alleles flipped (A1/A2 swapped)"
  } else {
    # alleles don't match (strand issue or ambiguous)
    ov$direction_note[i] <- sprintf("allele mismatch IA:%s/%s vs AAA:%s/%s", a1i,a2i,a1a,a2a)
    next
  }
  ov$direction_concordant[i] <- ifelse(conc, "concordant", "discordant")
}

# ------------------------------------------------------------
# Classify Type A/B/C
# ------------------------------------------------------------
ov[, overlap_type := NA_character_]
for (i in seq_len(nrow(ov))) {
  d <- ov$lead_distance[i]
  dir <- ov$direction_concordant[i]
  # Type A: lead distance < 200kb AND direction concordant
  # Type C: overlap but direction discordant
  # Type B: overlap but lead distance > 200kb (different causal variants)
  if (d <= 200000 && !is.na(dir) && dir == "concordant") {
    ov$overlap_type[i] <- "A"
  } else if (!is.na(dir) && dir == "discordant") {
    ov$overlap_type[i] <- "C"
  } else if (d > 200000) {
    ov$overlap_type[i] <- "B"
  } else if (d <= 200000 && is.na(dir)) {
    # close leads but unknown direction -> call A with note
    ov$overlap_type[i] <- "A?"
  } else {
    ov$overlap_type[i] <- "B"
  }
}

# sort by joint significance (product of P)
ov[, joint_p := ia_leadP * aaa_leadP]
setorder(ov, joint_p)

cat("\nOverlap summary by type:\n")
print(table(ov$overlap_type, useNA = "ifany"))
cat("\nDirection concordance:\n")
print(table(ov$direction_concordant, useNA = "ifany"))

# Output
setcolorder(ov, c("chr","ia_start","ia_end","aaa_start","aaa_end",
                  "ov_start","ov_end","ov_width","lead_distance",
                  "ia_leadSNP","ia_leadPos","ia_leadP","ia_leadBETA","ia_leadSE",
                  "aaa_leadSNP","aaa_leadPos","aaa_leadP","aaa_leadBETA","aaa_leadSE",
                  "direction_concordant","direction_note","overlap_type","joint_p"))
fwrite(ov, file.path(RESULTS_DIR, "locus_overlap.tsv"), sep = "\t")

cat("\nTop 15 overlapping loci:\n")
print(ov[1:min(15,.N), .(chr, lead_distance, ia_leadSNP, ia_leadP,
                          aaa_leadSNP, aaa_leadP, direction_concordant, overlap_type)])
cat("\nDone. Overlap file written.\n")
