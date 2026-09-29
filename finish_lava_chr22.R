# ============================================================
# Finish remaining loci 2401-2495 (chr22)
# ============================================================
library(LAVA)
library(data.table)

ref.prefix   <- "D:/IA_AAA_analysis/data/ld_reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.merged"
loc.file     <- "D:/IA_AAA_analysis/data/ld_reference/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile"
info.file    <- "D:/IA_AAA_analysis/data/lava_input.info.txt"
out.dir      <- "D:/IA_AAA_analysis/results/lava"

loci <- read.loci(loc.file)
cat("Total loci:", nrow(loci), "\n")

cat("Processing input...\n")
input <- process.input(
  input.info.file = info.file,
  sample.overlap.file = NULL,
  ref.prefix = ref.prefix,
  phenos = c("IA", "AAA", "SBP", "SmokingInit")
)

univ.thresh <- 1e-4

# Load existing results
existing_univ <- fread(file.path(out.dir, "lava_univariate_v2.tsv"))
existing_bivar <- fread(file.path(out.dir, "lava_bivar_raw_v2.tsv"))
existing_pcor <- fread(file.path(out.dir, "lava_pcor_conditional_v2.tsv"))
cat("Existing univ loci:", length(unique(existing_univ$locus_id)), "\n")
cat("Existing bivar rows:", nrow(existing_bivar), "\n")

# Process loci 2401-2495
start.idx <- 2401
end.idx <- 2495
cat("Processing loci", start.idx, "to", end.idx, "\n")

univ.list <- list()
bivar.list <- list()
pcor.list <- list()

get.p <- function(df, ph) {
  if (is.null(df)) return(NA_real_)
  v <- df$p[df$phen == ph]
  if (length(v) == 0) return(NA_real_)
  return(v[1])
}

for (i in start.idx:end.idx) {
  if (i %% 10 == 0) cat("..locus", i, "\n")

  locus <- tryCatch(
    process.locus(loci[i,], input),
    error = function(e) { cat("Error at locus", i, ":", e$message, "\n"); NULL }
  )
  if (is.null(locus)) next

  loc.info <- data.frame(
    locus_id = locus$id, chr = locus$chr, start = locus$start,
    stop = locus$stop, n_snps = locus$n.snps, n_pcs = locus$K
  )

  univ.df <- tryCatch(run.univ(locus), error = function(e) NULL)
  if (!is.null(univ.df)) univ.list[[length(univ.list)+1]] <- cbind(loc.info, univ.df)

  ia.p  <- get.p(univ.df, "IA")
  aaa.p <- get.p(univ.df, "AAA")

  if (!is.na(ia.p) && !is.na(aaa.p) && ia.p < univ.thresh && aaa.p < univ.thresh) {
    bv <- tryCatch(run.bivar(locus, phenos = c("IA", "AAA")), error = function(e) NULL)
    if (!is.null(bv)) bivar.list[[length(bivar.list)+1]] <- cbind(loc.info, bv)

    pc <- tryCatch(run.pcor(locus, target = c("IA", "AAA"),
                            phenos = c("SBP", "SmokingInit")), error = function(e) NULL)
    if (!is.null(pc)) pcor.list[[length(pcor.list)+1]] <- cbind(loc.info, pc)
  }
}

# Append to existing results
if (length(univ.list) > 0) {
  new_univ <- do.call(rbind, univ.list)
  all_univ <- rbind(existing_univ, new_univ)
  write.table(all_univ, file.path(out.dir, "lava_univariate_v2.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  cat("New univ rows:", nrow(new_univ), " Total:", nrow(all_univ), "\n")
}

if (length(bivar.list) > 0) {
  new_bivar <- do.call(rbind, bivar.list)
  all_bivar <- rbind(existing_bivar, new_bivar)
  write.table(all_bivar, file.path(out.dir, "lava_bivar_raw_v2.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  cat("New bivar rows:", nrow(new_bivar), " Total:", nrow(all_bivar), "\n")
}

if (length(pcor.list) > 0) {
  new_pcor <- do.call(rbind, pcor.list)
  all_pcor <- rbind(existing_pcor, new_pcor)
  write.table(all_pcor, file.path(out.dir, "lava_pcor_conditional_v2.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  cat("New pcor rows:", nrow(new_pcor), " Total:", nrow(all_pcor), "\n")
}

cat("\n=== DONE. Total univ loci:", length(unique(all_univ$locus_id)), "/ 2495 ===\n")
