# ============================================================
# Diagnose existing LAVA results before re-running
# ============================================================
library(data.table)

out.dir <- "D:/IA_AAA_analysis/results/lava"

# Load existing results
univ <- fread(file.path(out.dir, "lava_univariate.tsv"))
bivar <- fread(file.path(out.dir, "lava_bivar_raw.tsv"))
pcor <- fread(file.path(out.dir, "lava_pcor_conditional.tsv"))

cat("=== EXISTING RESULTS SUMMARY ===\n")
cat("Univariate: ", nrow(univ), " rows, ", length(unique(univ$locus_id)), " unique loci\n")
cat("Bivar raw:  ", nrow(bivar), " rows, ", length(unique(bivar$locus_id)), " unique loci\n")
cat("Pcor cond:  ", nrow(pcor), " rows, ", length(unique(pcor$locus_id)), " unique loci\n")

# Which loci are missing?
all.loci <- 1:2495
done.loci <- unique(univ$locus_id)
missing.loci <- setdiff(all.loci, done.loci)
cat("\nMissing loci: ", length(missing.loci), "\n")
cat("Missing locus range: ", min(missing.loci), " - ", max(missing.loci), "\n")

# Check chr3 rho = 1.00
cat("\n=== CHR3 RHO CHECK ===\n")
chr3.bivar <- bivar[chr == 3]
cat("Chr3 bivar loci: ", nrow(chr3.bivar), "\n")
rho1 <- chr3.bivar[rho >= 0.999]
cat("Chr3 with rho >= 0.999: ", nrow(rho1), "\n")
if (nrow(rho1) > 0) {
  print(rho1[, .(locus_id, chr, start, stop, rho, rho.lower, rho.upper, p, n_snps, n_pcs)])
}

# Check all rho near 1
cat("\nAll loci with rho >= 0.999:\n")
rho.all <- bivar[rho >= 0.999]
print(rho.all[, .(locus_id, chr, start, stop, rho, p)])

# Check NA in pcor
cat("\n=== PCOR NA DIAGNOSTICS ===\n")
cat("Total pcor rows: ", nrow(pcor), "\n")
cat("Rows with NA p: ", sum(is.na(pcor$p)), "\n")
cat("Rows with NA pcor: ", sum(is.na(pcor$pcor)), "\n")
cat("Rows with NA ci.lower: ", sum(is.na(pcor$ci.lower)), "\n")
cat("Rows with NA ci.upper: ", sum(is.na(pcor$ci.upper)), "\n")

# Look at some NA examples
na.rows <- pcor[is.na(p)]
cat("\nNA p rows - first 10:\n")
if (nrow(na.rows) > 0) {
  print(head(na.rows, 10))
  cat("\nNA rows by chr:\n")
  print(na.rows[, .N, by = chr][order(chr)])
}

# Check the non-NA pcor rows
ok.rows <- pcor[!is.na(p)]
cat("\nNon-NA pcor rows: ", nrow(ok.rows), "\n")
cat("Significant at Bonferroni (0.05/2495): ", sum(ok.rows$p < 0.05/2495), "\n")
cat("Significant at FDR (0.05): ", sum(p.adjust(ok.rows$p, method="fdr") < 0.05), "\n")

# Check bivar significance
cat("\n=== BIVAR SIGNIFICANCE (old, p<0.05 threshold) ===\n")
cat("Bivar rows: ", nrow(bivar), "\n")
cat("Significant at Bonferroni (0.05/2495): ", sum(bivar$p < 0.05/2495, na.rm=TRUE), "\n")
cat("Significant at FDR (0.05): ", sum(p.adjust(bivar$p, method="fdr") < 0.05, na.rm=TRUE), "\n")

# Check univariate p-value distribution
cat("\n=== UNIVARIATE P-VALUE DISTRIBUTION ===\n")
ia.univ <- univ[phen == "IA"]
aaa.univ <- univ[phen == "AAA"]
cat("IA: ", nrow(ia.univ), " loci with h2 estimate\n")
cat("IA p < 1e-4: ", sum(ia.univ$p < 1e-4, na.rm=TRUE), "\n")
cat("IA p < 0.05: ", sum(ia.univ$p < 0.05, na.rm=TRUE), "\n")
cat("AAA: ", nrow(aaa.univ), " loci with h2 estimate\n")
cat("AAA p < 1e-4: ", sum(aaa.univ$p < 1e-4, na.rm=TRUE), "\n")
cat("AAA p < 0.05: ", sum(aaa.univ$p < 0.05, na.rm=TRUE), "\n")

# How many loci pass BOTH at p<1e-4?
ia.sig <- ia.univ[p < 1e-4, locus_id]
aaa.sig <- aaa.univ[p < 1e-4, locus_id]
both.sig <- intersect(ia.sig, aaa.sig)
cat("\nLoci passing BOTH IA and AAA at p<1e-4: ", length(both.sig), "\n")
cat("(Previously at p<0.05: 748 bivar rows)\n")

cat("\n=== DONE DIAGNOSTICS ===\n")
