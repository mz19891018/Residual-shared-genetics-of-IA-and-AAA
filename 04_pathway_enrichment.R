# ============================================================
# 04_pathway_enrichment.R
# Pathway enrichment for IA / AAA / overlap locus gene sets
# Uses Fisher's exact test with custom vascular gene sets + GO from biomaRt
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

RESULTS_DIR <- "D:/IA_AAA_analysis/results"

# Load genes
genes <- fread(file.path(RESULTS_DIR, "grch37_protein_coding_genes.tsv"))
setnames(genes, "gene", "hgnc")
universe <- unique(genes$hgnc)
cat("Gene universe (protein-coding):", length(universe), "\n")

# Load loci
ia_loci <- fread(file.path(RESULTS_DIR, "ia_significant_loci.tsv"))
aaa_loci <- fread(file.path(RESULTS_DIR, "aaa_significant_loci.tsv"))
ov <- fread(file.path(RESULTS_DIR, "locus_overlap.tsv"))

# ------------------------------------------------------------
# Collect genes within locus boundaries
# ------------------------------------------------------------
locus_genes <- function(loci_dt, gene_dt) {
  all_genes <- c()
  for (i in seq_len(nrow(loci_dt))) {
    r <- loci_dt[i]
    g <- gene_dt[chr == as.character(r$chr) & end >= r$start & start <= r$end, hgnc]
    all_genes <- c(all_genes, g)
  }
  unique(all_genes)
}

ia_genes <- locus_genes(ia_loci, genes)
aaa_genes <- locus_genes(aaa_loci, genes)
# overlap loci = union of IA and AAA genes in overlapping regions
ov_ia_loci <- ia_loci[paste(chr,start,end) %in% paste(ov$chr,ov$ia_start,ov$ia_end)]
ov_aaa_loci <- aaa_loci[paste(chr,start,end) %in% paste(ov$chr,ov$aaa_start,ov$aaa_end)]
ov_genes <- unique(c(locus_genes(ov_ia_loci, genes), locus_genes(ov_aaa_loci, genes)))
# concordant direction overlap genes (Type A)
conc_idx <- which(ov$direction_concordant == "concordant")
conc_ia <- ia_loci[paste(chr,start,end) %in%
                   paste(ov$chr[conc_idx], ov$ia_start[conc_idx], ov$ia_end[conc_idx])]
conc_aaa <- aaa_loci[paste(chr,start,end) %in%
                     paste(ov$chr[conc_idx], ov$aaa_start[conc_idx], ov$aaa_end[conc_idx])]
conc_genes <- unique(c(locus_genes(conc_ia, genes), locus_genes(conc_aaa, genes)))

cat("IA locus genes:", length(ia_genes), "\n")
cat("AAA locus genes:", length(aaa_genes), "\n")
cat("Overlap locus genes:", length(ov_genes), "\n")
cat("Concordant overlap genes:", length(conc_genes), "\n")

# ------------------------------------------------------------
# Custom vascular gene sets (curated from known biology)
# ------------------------------------------------------------
gene_sets <- list(
  "ECM_organization" = c("COL1A1","COL1A2","COL3A1","COL4A1","COL4A2","COL5A1","COL5A2",
    "COL6A1","COL6A2","COL6A3","COL18A1","ELN","FBN1","FBN2","LTBP1","LTBP2","LTBP3",
    "LTBP4","EMILIN1","EMILIN2","MFAP5","VCAN","ACAN","BGN","DCN","LUM","FN1",
    "LAMA1","LAMA2","LAMA3","LAMA4","LAMA5","LAMB1","LAMB2","LAMC1","HSPG2",
    "MATN1","MATN2","MATN3","THBS1","THBS2","THBS3","THBS4","COMP","SPARC","SPARCL1",
    "TNC","TNXB","FBLN1","FBLN2","FBLN5","EFEMP1","EFEMP2","LOX","LOXL1","LOXL2",
    "LOXL3","LOXL4","PLOD1","PLOD2","PLOD3","MMP2","MMP9","TIMP1","TIMP2","TIMP3"),

  "VSMC_contraction" = c("ACTA2","MYH11","MYLK","CNN1","CNN2","TAGLN","CALD1","SMTN",
    "TPM1","TPM2","ACTG2","MYL9","MYL6","PPP1R12A","ROCK1","ROCK2","MYOCD","SRF","ELN",
    "ACTN1","ACTN4","FLNA","FLNB","MYH9","MYL6B"),

  "Endothelial_function" = c("CDH5","PECAM1","VWF","ENG","KDR","FLT1","TEK","TIE1","ESAM",
    "CLDN5","OCLN","TJP1","TJP2","JAM2","JAM3","CDH2","VEZF1","SOX18","ERG","FLI1",
    "SOX17","NOS3","EDN1","EDNRA","EDNRB","ACE","ACE2","AGT","NPPB","VEGFA","VEGFB",
    "VEGFC","FLT4","NRP1","NRP2","ANGPT1","ANGPT2","EPHB4","EFNB2","DLL4","NOTCH4",
    "JAG1","HEY1","HEY2","ETV2","ETS1"),

  "TGFbeta_signaling" = c("TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","TGFBR3","SMAD2",
    "SMAD3","SMAD4","SMAD7","BMP2","BMP4","BMP6","BMP7","BMPR1A","BMPR2","ACVRL1",
    "ENG","LEFTY1","LEFTY2","NOG","CHRD","FST","LTBP1","LTBP2","LTBP4","SERPINE1","CTGF",
    "BMP10","GDF2","GDF5","MYOCD","SRF"),

  "Inflammation" = c("TNF","IL1B","IL6","IL6R","IL10","NFKB1","RELA","RELB","NFKBIA",
    "TLR2","TLR4","CD14","MYD88","NLRP3","IL18","CCL2","CCL5","CXCL8","CXCL10",
    "ICAM1","VCAM1","SELE","SELP","CRP","PTGS2","MMP1","MMP3","MMP12","IL1A","IL1RN",
    "IFNG","STAT1","STAT3","IRF1","IRF8"),

  "Lipid_metabolism" = c("LDLR","PCSK9","SORT1","APOB","APOE","APOA1","APOA2","APOA4",
    "APOC1","APOC2","APOC3","LPL","LIPC","LIPG","CETP","LCAT","ABCA1","ABCG1","ABCG5",
    "ABCG8","HMGCR","HMGCS1","NPC1L1","FADS1","FADS2","LIPA","SOAT1","APOE","LDLRAP1",
    "PCSK9","ANGPTL3","ANGPTL4","LPL","PLTP","SCARB1","CD36","VLDLR","LRP1","LRP2"),

  "Monogenic_aortopathy" = c("FBN1","TGFBR1","TGFBR2","SMAD3","COL3A1","ACTA2","MYH11",
    "MYLK","PRKG1","FLNA","LOX","MAT2A","FOXE3","NOTCH1","ELN","SLC2A10","COL5A1",
    "COL5A2","SMAD2","TGFB2","TGFB3","BGN","MFAP5","SLC2A10","TREX1","COL4A1","COL4A2"),

  "Coagulation_hemostasis" = c("F2","F3","F5","F7","F8","F9","F10","F11","F12","F13A1",
    "F13B","VWF","FGA","FGB","FGG","SERPINC1","SERPIND1","PROC","PROS1","PROZ","THBD",
    "PLAT","PLAU","SERPINE1","PLAT","ITGA2B","ITGB3","GP1BA","GP1BB","GP6","P2RY12",
    "PTGS1","NOS3","NOS2"),

  "Extracellular_matrix_mMP" = c("MMP1","MMP2","MMP3","MMP7","MMP8","MMP9","MMP10","MMP12",
    "MMP13","MMP14","MMP15","MMP16","MMP19","MMP20","MMP24","TIMP1","TIMP2","TIMP3",
    "TIMP4","ADAMTS1","ADAMTS4","ADAMTS5","ADAMTS13","PLAT","PLAU","SERPINE1"),

  "Endothelial_NO_signaling" = c("NOS3","NOS2","NOS1","AGT","ACE","ACE2","AGTR1","AGTR2",
    "REN","EDN1","EDNRA","EDNRB","NPPA","NPPB","NPPC","NPRA","NPRB","EDN1","EDN2","EDN3",
    "VEGFA","KDR","NOS3","CALCRL","RAMP1","RAMP2","RAMP3","ADM","ADM2")
)

# Filter gene sets to universe
gene_sets <- lapply(gene_sets, function(gs) intersect(gs, universe))

# ------------------------------------------------------------
# Fisher's exact test
# ------------------------------------------------------------
fisher_enrich <- function(query_genes, gene_sets, universe) {
  query <- intersect(query_genes, universe)
  n_query <- length(query)
  n_univ <- length(universe)
  res <- list()
  for (nm in names(gene_sets)) {
    gs <- gene_sets[[nm]]
    a <- length(intersect(query, gs))          # in both
    b <- length(setdiff(gs, query))             # in set, not query
    c <- n_query - a                            # in query, not set
    d <- n_univ - a - b - c                     # neither
    if (a == 0) {
      p <- 1
      or <- 0
    } else {
      ft <- fisher.test(matrix(c(a,b,c,d), nrow=2), alternative = "greater")
      p <- ft$p.value
      or <- unname(ft$estimate)
    }
    # expected overlap
    expected <- n_query * length(gs) / n_univ
    res[[nm]] <- data.table(
      pathway = nm,
      n_in_pathway = length(gs),
      n_in_query = n_query,
      overlap_count = a,
      expected = round(expected, 2),
      enrichment_ratio = round(a / max(expected, 0.01), 2),
      odds_ratio = round(or, 3),
      p_value = signif(p, 3),
      overlap_genes = paste(intersect(query, gs), collapse = ",")
    )
  }
  rbindlist(res)
}

cat("\n=== Enrichment: IA all locus genes ===\n")
ia_enr <- fisher_enrich(ia_genes, gene_sets, universe)
ia_enr[, gene_set := "IA_all_loci"]
print(ia_enr[order(p_value), .(pathway, overlap_count, enrichment_ratio, p_value)])

cat("\n=== Enrichment: AAA all locus genes ===\n")
aaa_enr <- fisher_enrich(aaa_genes, gene_sets, universe)
aaa_enr[, gene_set := "AAA_all_loci"]
print(aaa_enr[order(p_value), .(pathway, overlap_count, enrichment_ratio, p_value)])

cat("\n=== Enrichment: Overlap locus genes ===\n")
ov_enr <- fisher_enrich(ov_genes, gene_sets, universe)
ov_enr[, gene_set := "IA_AAA_overlap"]
print(ov_enr[order(p_value), .(pathway, overlap_count, enrichment_ratio, p_value)])

cat("\n=== Enrichment: Concordant direction overlap genes ===\n")
conc_enr <- fisher_enrich(conc_genes, gene_sets, universe)
conc_enr[, gene_set := "Overlap_concordant"]
print(conc_enr[order(p_value), .(pathway, overlap_count, enrichment_ratio, p_value)])

# Combine and adjust p-values
all_enr <- rbind(ia_enr, aaa_enr, ov_enr, conc_enr)
all_enr[, FDR := p.adjust(p_value, method = "BH"), by = gene_set]
setcolorder(all_enr, c("gene_set","pathway","n_in_pathway","n_in_query",
                        "overlap_count","expected","enrichment_ratio",
                        "odds_ratio","p_value","FDR","overlap_genes"))
setorder(all_enr, gene_set, p_value)

fwrite(all_enr, file.path(RESULTS_DIR, "pathway_enrichment.tsv"), sep = "\t")

# Save gene sets for report
cat("\nOverlap locus genes:", paste(ov_genes, collapse=", "), "\n")
cat("\nConcordant overlap genes:", paste(conc_genes, collapse=", "), "\n")
cat("\nPathway enrichment saved.\n")
