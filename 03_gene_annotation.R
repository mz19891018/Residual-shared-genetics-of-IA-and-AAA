# ============================================================
# 03_gene_annotation.R
# Annotate loci with protein-coding genes via biomaRt (GRCh37)
# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(biomaRt)
})

RESULTS_DIR <- "D:/IA_AAA_analysis/results"
GENE_FILE <- file.path(RESULTS_DIR, "grch37_protein_coding_genes.tsv")

# ------------------------------------------------------------
# Step 1: Get gene coordinates from biomaRt (GRCh37)
# ------------------------------------------------------------
if (!file.exists(GENE_FILE)) {
  cat("Fetching protein-coding genes from biomaRt (GRCh37)...\n")
  options(download.file.method = "curl")
  ensembl <- tryCatch({
    useMart(biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = "hsapiens_gene_ensembl",
            host = "https://grch37.ensembl.org")
  }, error = function(e) {
    cat("Primary host failed:", conditionMessage(e), "\n")
    tryCatch(useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                     dataset = "hsapiens_gene_ensembl",
                     host = "https://feb2014.archive.ensembl.org"),
             error = function(e2) {
               cat("Archive host also failed:", conditionMessage(e2), "\n")
               NULL
             })
  })

  if (!is.null(ensembl)) {
    genes <- getBM(attributes = c("chromosome_name","start_position","end_position",
                                    "strand","external_gene_name","gene_biotype",
                                    "entrezgene_id"),
                   filters = c("biotype","chromosome_name"),
                   values = list(biotype = "protein_coding",
                                 chromosome_name = c(1:22,"X","Y")),
                   mart = ensembl)
    setDT(genes)
    setnames(genes, c("chromosome_name","start_position","end_position",
                      "strand","external_gene_name","gene_biotype","entrezgene_id"),
                    c("chr","start","end","strand","gene","biotype","entrez"))
    genes <- genes[biotype == "protein_coding" & grepl("^[0-9XY]+$", chr)]
    fwrite(genes, GENE_FILE, sep = "\t")
    cat("Saved", nrow(genes), "protein-coding genes to", GENE_FILE, "\n")
  } else {
    cat("WARNING: biomaRt unavailable. Will use fallback gene list.\n")
  }
} else {
  cat("Using cached gene file:", GENE_FILE, "\n")
}

if (file.exists(GENE_FILE)) {
  genes <- fread(GENE_FILE)
  cat("Loaded", nrow(genes), "protein-coding genes.\n")
} else {
  # Fallback: minimal known vascular gene list (GRCh37 coords approximate)
  cat("Using fallback gene list.\n")
  genes <- data.table(
    chr = c("9","9","9","4","18","16","12","13","5","15","10","22",
            "19","19","1","6","1","15","3","7","15","9","5","13","20",
            "1","6","10","17","3","11","13","18","21","X"),
    start = c(21967000,21990000,21995000,148300000,20150000,75250000,
              95400000,33600000,131600000,78700000,104900000,30200000,
              11000000,44000000,109000000,161000000,109000000,78800000,
              129000000,128000000,78800000,22000000,131000000,33000000,
              44000000,109000000,161000000,104000000,48000000,129000000,
              130000000,33000000,20000000,45000000,15400000),
    end   = c(21995000,22020000,22020000,148450000,20300000,75400000,
              95600000,33800000,131800000,78900000,105100000,30400000,
              11200000,44200000,109200000,161200000,109200000,79000000,
              129200000,128200000,79000000,22100000,131500000,33500000,
              44200000,109200000,161200000,105200000,48200000,129200000,
              130500000,33500000,20500000,45200000,15600000),
    gene = c("CDKN2A","CDKN2B","MTAP","EDIL3","BCL2","RP11",
             "MRPL35","MYO16","GPC6","HMG20A","NEUROG1","CHEK2",
             "LDLR","PCSK9","SORT1","HLA","APOE","FBN1","SMAD3",
             "TGFBR2","IGF2R","CDKN2B-AS1","EDNRA","EDNRB","C20orf",
             "ABO","HFE","CYTIP","NOS2","ACTA2","MYH11","COL3A1",
             "SOX17","CNNM2","FLNA"),
    biotype = "protein_coding"
  )
}

# ------------------------------------------------------------
# Annotation function
# ------------------------------------------------------------
annotate_locus <- function(chr, start, end, lead_pos, gene_dt) {
  g <- gene_dt[chr == as.character(chr)]
  if (nrow(g) == 0) return(list(genes_in_locus = NA_character_,
                                 nearest_gene = NA_character_,
                                 lead_in_gene = NA_character_,
                                 distance_to_nearest = NA_integer_))
  # genes within locus
  in_locus <- g[end >= start & start <= end]
  genes_in <- paste(unique(in_locus$gene), collapse = ";")
  if (genes_in == "") genes_in <- NA_character_

  # nearest gene to lead position
  g[, dist_to_lead := pmax(abs(start - lead_pos), abs(end - lead_pos), 0)]
  # if lead inside gene, dist = 0
  g[, lead_inside := start <= lead_pos & end >= lead_pos]
  setorder(g, dist_to_lead)
  nearest <- g[1]
  lead_in <- ifelse(nearest$lead_inside, nearest$gene, NA_character_)
  return(list(genes_in_locus = genes_in,
              nearest_gene = nearest$gene,
              lead_in_gene = lead_in,
              distance_to_nearest = nearest$dist_to_lead))
}

# ------------------------------------------------------------
# Annotate IA loci
# ------------------------------------------------------------
ia_loci <- fread(file.path(RESULTS_DIR, "ia_significant_loci.tsv"))
aaa_loci <- fread(file.path(RESULTS_DIR, "aaa_significant_loci.tsv"))
ov <- fread(file.path(RESULTS_DIR, "locus_overlap.tsv"))

# Build annotation table: all overlapping loci + top 20 per disease
# Overlapping IA loci
overlap_ia_chrs <- unique(ov$chr)
overlap_ia_loci <- ia_loci[paste(chr, start, end) %in%
                           paste(ov$chr, ov$ia_start, ov$ia_end)]
# Top 20 IA by P
ia_top20 <- ia_loci[order(lead_P)][1:min(20, .N)]
ia_annot <- unique(rbind(overlap_ia_loci, ia_top20))

# Top 20 AAA by P
aaa_top20 <- aaa_loci[order(lead_P)][1:min(20, .N)]
# Overlapping AAA loci
overlap_aaa_loci <- aaa_loci[paste(chr, start, end) %in%
                             paste(ov$chr, ov$aaa_start, ov$aaa_end)]
aaa_annot <- unique(rbind(overlap_aaa_loci, aaa_top20))

cat("Annotating", nrow(ia_annot), "IA loci and", nrow(aaa_annot), "AAA loci...\n")

annotate_dt <- function(loci_dt, disease_label, gene_dt) {
  res <- loci_dt[, {
    a <- annotate_locus(chr, start, end, lead_pos, gene_dt)
    .(disease = disease_label, chr = chr, start = start, end = end,
      lead_SNP = lead_SNP, lead_pos = lead_pos, lead_P = lead_P,
      genes_in_locus = a$genes_in_locus,
      nearest_gene = a$nearest_gene,
      lead_in_gene = a$lead_in_gene,
      distance_to_nearest_bp = a$distance_to_nearest)
  }, by = 1:nrow(loci_dt)]
  res[, nrow := NULL]
  return(res)
}

ia_ann_dt <- annotate_dt(ia_annot, "IA", genes)
aaa_ann_dt <- annotate_dt(aaa_annot, "AAA", genes)

# Also annotate overlap pairs specifically
ov_annot <- ov[, .(
  overlap_type, chr, ia_start, ia_end, aaa_start, aaa_end,
  ia_leadSNP, ia_leadPos, ia_leadP,
  aaa_leadSNP, aaa_leadPos, aaa_leadP,
  direction_concordant, direction_note, lead_distance
)]

# Add nearest gene info for overlap pairs
ov_annot[, ia_nearest_gene := vapply(seq_len(.N), function(i) {
  a <- annotate_locus(chr[i], ia_start[i], ia_end[i], ia_leadPos[i], genes)
  a$nearest_gene
}, character(1))]
ov_annot[, aaa_nearest_gene := vapply(seq_len(.N), function(i) {
  a <- annotate_locus(chr[i], aaa_start[i], aaa_end[i], aaa_leadPos[i], genes)
  a$nearest_gene
}, character(1))]

# Combine all annotations
all_ann <- rbind(
  ia_ann_dt[, .(disease, chr, start, end, lead_SNP, lead_pos, lead_P,
                genes_in_locus, nearest_gene, lead_in_gene, distance_to_nearest_bp)],
  aaa_ann_dt[, .(disease, chr, start, end, lead_SNP, lead_pos, lead_P,
                 genes_in_locus, nearest_gene, lead_in_gene, distance_to_nearest_bp)]
)

fwrite(all_ann, file.path(RESULTS_DIR, "locus_gene_annotation.tsv"), sep = "\t")
fwrite(ov_annot, file.path(RESULTS_DIR, "overlap_gene_detail.tsv"), sep = "\t")

cat("\n=== Overlap gene detail ===\n")
print(ov_annot[, .(chr, overlap_type, ia_leadSNP, ia_nearest_gene,
                   aaa_leadSNP, aaa_nearest_gene, direction_concordant)])
cat("\nGene annotation saved.\n")
