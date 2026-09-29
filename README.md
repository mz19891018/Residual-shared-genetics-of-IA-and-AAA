# Residual-shared-genetics-of-IA-and-AAA
The shared loci at 9p21.3 and 18q11 likely act through non-coding mechanisms not captured by bulk arterial gene expression.
IA–AAA Shared Genetic Architecture: Analysis Code README
Project Overview
This repository contains the analysis code for the study:
"Shared genetic architecture between intracranial aneurysm and abdominal aortic aneurysm"
All data are stored in D:\IA_AAA_analysis\data\ (D drive), R is installed on C drive (C:\Program Files\R\R-4.6.0\bin\Rscript.exe), and results are written to D:\IA_AAA_analysis\results\.

Directory Structure
D:\IA_AAA_analysis\
├── data\                          # Raw and processed input data
│   ├── gwas\
│   │   ├── ia_2020\               # Bakker 2020 IA GWAS (Stage 1)
│   │   ├── aaa_2023\              # AAAGen 2023 AAA GWAS
│   │   ├── sbp\                   # ICBP 2011 SBP
│   │   ├── smoking\               # GSCAN smoking initiation
│   │   ├── cad\                   # CARDIoGRAMplusC4D CAD
│   │   ├── lipids\                # GLGC lipids
│   │   └── height\                # GIANT height
│   ├── ldsc\                      # LDSC LD scores (UKB baseline v2.2)
│   ├── gtex\                      # GTEx v8 EUR eQTL
│   └── singlecell\                # GSE166676, GSE193533 RAW tar files
├── tools\                         # All analysis scripts (this directory)
├── results\                       # Output tables, figures, and intermediate files
│   ├── ldsc\
│   ├── genomic_sem\
│   ├── lava\
│   ├── coloc\
│   ├── susie\
│   ├── twas\
│   ├── mr\
│   ├── singlecell\
│   └── *.png, *.docx              # Final figures and reports
└── Manuscript_v6_with_supplements.docx


Analysis Pipeline (Run Order)
Step 1: Data Preparation & Quality Control
Script
Language
Purpose
qc_munge.py
Python
Munge GWAS summary stats for LDSC (align alleles, filter MAF/INFO)
fix_alleles.py
Python
Fix strand flips and allele mismatches between GWAS
extract_ldsc_formals.R
R
Extract LDSC input from sumstats

Key inputs: GWAS sumstats in data/gwas/
Key outputs: LDSC-ready .sumstats.gz files in results/ldsc/

Step 2: LDSC Genetic Correlation
Script
Language
Purpose
ldsc.py
Python
Official LDSC (Bulik-Sullivan) — rg and h²
munge_sumstats.py
Python
Official LDSC munging
ldscore.py
Python
Official LDSC LD score calculation

Run:
python tools/ldsc.py --rg ia_2020.sumstats.gz,aaa_2023.sumstats.gz --ref-ld-chr data/ldsc/ --w-ld-chr data/ldsc/

Key outputs: results/ldsc/ia_aaa_rg.log → rg=0.349 (SE 0.064, P=5.0e-8)

Step 3: Genomic SEM (Conditional Genetic Correlation)
Script
Language
Purpose
run_genomic_sem.R
R
Base model: IA-AAA rg conditioning on SBP+smoking
run_genomic_sem_extended.R
R
Extended models: +LDL (Model 5), +CAD (Model 6)
plot_forest_extended.R
R
Forest plot of residual rg across 6 models

Key outputs:
results/genomic_sem/extended_conditional_rg_forest.png (Figure 2)
Residual rg: Raw=0.331, +SBP=0.330, +Smoking=0.269, +SBP+Smoking=0.260, +LDL=0.266, +CAD=0.254

Step 4: LAVA Local Genetic Correlation
Script
Language
Purpose
run_lava_v2.R
R
LAVA bivariate local rg across 2,495 LD blocks
postprocess_lava_v2.R
R
Filter to 69 testable blocks, apply Bonferroni (0.05/69)
summarize_lava.R
R
Summarize significant loci
diagnose_lava_v2.R
R
Boundary flag diagnostics

Key outputs:
results/lava/lava_bivar_raw_v2_clean.tsv (69 blocks)
Two significant loci: 9p21.3 (ρ=0.727, P=4.0e-13), 18q11.2 (ρ=0.675, P=4.9e-7)
results/Fig3_LAVA_manhattan.png (Figure 3)

Step 5: Locus Overlap & Fine-Mapping
Script
Language
Purpose
01_significant_loci.R
R
Extract genome-wide significant loci for IA and AAA
02_locus_overlap.R
R
Identify overlapping loci within 500 kb
run_susie_coloc.R
R
SuSiE-coloc fine-mapping at overlapping loci
plot_susie.R
R
PIP plots for SuSiE results

Key outputs: results/susie/coloc_susie_ia_aaa_summary.tsv

Step 6: eQTL Colocalization
Script
Language
Purpose
extract_mtap_eqtl.py
Python
Extract MTAP eQTL SNPs from GTEx v8 Artery Aorta
run_mtap_coloc.R
R
coloc.abf for MTAP eQTL vs IA/AAA GWAS
run_cdkn2a_eqtl_coloc.R
R
coloc.abf for CDKN2A eQTL vs IA/AAA GWAS
prepare_coloc_data.py
Python
Prepare aligned SNP sets for coloc

Key outputs:
results/coloc/MTAP_coloc_results.tsv (PP.H4<0.03 — no colocalization)
results/Fig4C_9p21_locus_MTAP_eQTL.png (Figure 4C)

Step 7: TWAS (Transcriptome-wide Association)
Script
Language
Purpose
run_twas.R
R
FUSION TWAS using GTEx v8 artery weights
twas_fast.R
R
Rapid scan across all genes
twas_postprocess.R
R
Identify shared genes between IA and AAA
FUSION.assoc_test.R
R
Official FUSION association test
FUSION.compute_weights.R
R
Official FUSION weight computation

Key outputs:
results/twas/twas_shared_genes_mapped.tsv (20 shared genes)
results/twas/twas_volcano_*.png (Supplementary Figure S3)

Step 8: Mendelian Randomization
Script
Language
Purpose
mr_analysis.R
R
Two-sample MR (IVW, MR-Egger, weighted median)
run_mr_raps.py
Python
MR-RAPS robust adjustment for SBP→AAA
run_mr_raps.R
R
MR-RAPS R implementation
plot_mr_raps.py
Python
Forest plots with RAPS results

Key outputs:
results/mr/mr_summary.tsv (Supplementary Table S1)
results/mr/mr_scatter_*.png (Supplementary Figures S1–S2)

Step 9: Single-Cell Analysis
Script
Language
Purpose
run_human_seurat.R
R
Human AAA Seurat clustering (GSE166676, 12,713 cells)
regenerate_umap.R
R
Regenerate UMAP with refined annotations
run_mouse_seurat.R
R
Mouse IA Seurat (GSE193533)
extract_cdkn2a.R
R
Extract CDKN2A expression by cell type from Seurat object
make_fig5_final.py
Python
Figure 5: UMAP + CDKN2A dot plots (human + mouse)

Key outputs:
results/singlecell/human_aaa_seurat_regenerated.rds
results/singlecell/cdkn2a_by_celltype_new.tsv
results/Fig5_human_UMAP_with_CDKN2A.png (Figure 5)

Step 10: Figure Generation
Script
Language
Purpose
redraw_figures.py
Python
Generate Figures 1, 3, 4C, 5
redraw_fig4c_final.py
Python
Figure 4C: 9p21.3 locus with real SNP scatter
redraw_fig2.py
Python
Figure 2: Genomic SEM forest (no internal legend)
make_fig5_final.py
Python
Figure 5: 3-panel horizontal layout

Final figures in results/:
Fig1_LDSC_heatmap.png — 9-trait LDSC correlation heatmap
genomic_sem/extended_conditional_rg_forest.png — Figure 2
Fig3_LAVA_manhattan.png — LAVA Manhattan with 0.05/69 threshold
Fig4C_9p21_locus_MTAP_eQTL.png — Real SNP locus plot + MTAP eQTL
Fig5_human_UMAP_with_CDKN2A.png — UMAP + CDKN2A expression dot plots

Step 11: Report & Manuscript Assembly
Script
Language
Purpose
generate_report_v7.py
Python
Generate analysis report Word document
insert_figures_final.py
Python
Insert figures into report
insert_supplements.py
Python
Insert supplementary tables and figures into manuscript


Software Requirements
R 4.6.0 (C:\Program Files\R\R-4.6.0\bin\Rscript.exe)
Seurat, LAVA, coloc, SusieR, TwoSampleMR, GenomicSEM, FUSION
Python 3.11+
pandas, numpy, matplotlib, scipy, python-docx
PLINK 2 (D:\IA_AAA_analysis\tools\plink2\plink2.exe)
LDSC (bundled in tools/)

Key Results Summary
Analysis
Key Finding
LDSC
IA–AAA rg = 0.349 (P=5.0e-8)
Genomic SEM
Residual rg after +CAD = 0.254 (not abolished)
LAVA
2 significant loci: 9p21.3 (ρ=0.727) and 18q11.2 (ρ=0.675)
Colocalization
MTAP eQTL does NOT colocalize (PP.H4<0.03)
TWAS
20 nominally shared genes (exploratory, in supplement)
MR
Smoking→IA β=0.63, →AAA β=0.52; SBP→IA β=0.059
scRNA-seq
12,713 human AAA cells; CDKN2A highest in VSMC


Notes
All numbers in the manuscript are traceable to result files in results/.
GWAS data are European ancestry only.
dbGaP-controlled data (phs001672) could not be accessed due to NIH policy; FinnGen and CVDKP summary statistics used instead.
Supplementary tables S1–S3 and figures S1–S3 are embedded in the final manuscript.
