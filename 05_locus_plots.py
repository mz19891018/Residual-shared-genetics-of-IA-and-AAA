# ============================================================
# 05_locus_plots.py
# LocusZoom-style plots for top 3 IA-AAA overlapping loci
# ============================================================
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os, gzip

RESULTS_DIR = r"D:\IA_AAA_analysis\results"
IA_FILE = r"D:\IA_AAA_analysis\data\gwas\ia_2020\IA.GWAS.BakkerMK.2020.sumstats.Stage_1.txt.gz"
AAA_FILE = r"D:\IA_AAA_analysis\data\gwas\aaa_2023\meta-AAAgen-final-sumstat.txt.gz"
GENE_FILE = os.path.join(RESULTS_DIR, "grch37_protein_coding_genes.tsv")

# Top 3 overlap loci by joint P (from locus_overlap.tsv)
# 1: chr9:22.1Mb  2: chr18:20.2Mb  3: chr15:78.8Mb
TOP_LOCI = [
    {"chr": 9,  "center": 22100000, "label": "chr9_22100k_CDKN2AB"},
    {"chr": 18, "center": 20250000, "label": "chr18_20250k_RBBP8"},
    {"chr": 15, "center": 78850000, "label": "chr15_78850k_CHRNA5"},
]
HALF_WINDOW = 500000  # +/- 500kb

# Load gene positions
genes = pd.read_csv(GENE_FILE, sep="\t")
genes["chr"] = genes["chr"].astype(str)

# ------------------------------------------------------------
# Read IA data (full, ~4.5M SNPs, manageable)
# ------------------------------------------------------------
print("Reading IA data...")
ia = pd.read_csv(IA_FILE, sep="\t",
                 usecols=["CHR","BP","SNP","BETA","SE","P"])
ia["CHR"] = ia["CHR"].astype(str)
print(f"  IA SNPs: {len(ia)}")

# ------------------------------------------------------------
# Read AAA data in chunks, extract needed regions
# ------------------------------------------------------------
print("Reading AAA data (chunked)...")
aaa_chunks = []
usecols = ["chr","pos","Effectsize","Effectsize_SD","pvalue"]
for chunk in pd.read_csv(AAA_FILE, sep="\t", usecols=usecols, chunksize=500000):
    # keep SNPs in any of our 3 regions
    for loc in TOP_LOCI:
        c = loc["chr"]; lo = loc["center"]-HALF_WINDOW; hi = loc["center"]+HALF_WINDOW
        sub = chunk[(chunk["chr"]==c) & (chunk["pos"]>=lo) & (chunk["pos"]<=hi)]
        if len(sub) > 0:
            aaa_chunks.append(sub)
aaa = pd.concat(aaa_chunks, ignore_index=True)
aaa["chr"] = aaa["chr"].astype(str)
aaa = aaa.rename(columns={"Effectsize":"BETA","Effectsize_SD":"SE","pvalue":"P"})
print(f"  AAA SNPs in target regions: {len(aaa)}")

# ------------------------------------------------------------
# Plotting function
# ------------------------------------------------------------
def plot_locus(loc):
    c = str(loc["chr"])
    lo = loc["center"] - HALF_WINDOW
    hi = loc["center"] + HALF_WINDOW

    ia_sub = ia[(ia["CHR"]==c) & (ia["BP"]>=lo) & (ia["BP"]<=hi)].copy()
    aaa_sub = aaa[(aaa["chr"]==c) & (aaa["pos"]>=lo) & (aaa["pos"]<=hi)].copy()

    ia_sub["neg_log10P"] = -np.log10(ia_sub["P"].clip(lower=1e-300))
    aaa_sub["neg_log10P"] = -np.log10(aaa_sub["P"].clip(lower=1e-300))

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True,
                              gridspec_kw={"height_ratios":[1,1]})

    # IA panel
    ax1 = axes[0]
    ax1.scatter(ia_sub["BP"]/1e6, ia_sub["neg_log10P"],
                s=8, alpha=0.5, c="steelblue", label="IA (Bakker 2020)")
    # highlight lead
    ia_lead_p = ia_sub["P"].min()
    ia_lead_pos = ia_sub.loc[ia_sub["P"].idxmin(), "BP"]
    ax1.scatter([ia_lead_pos/1e6], [-np.log10(ia_lead_p)],
                s=80, c="darkblue", marker="*", zorder=5, edgecolors="white")
    ax1.axhline(-np.log10(5e-8), color="red", ls="--", lw=0.8, label="P=5e-8")
    ax1.set_ylabel("-log10(P)  [IA]")
    ax1.set_title(f"Chromosome {c}: {lo/1e6:.1f}-{hi/1e6:.1f} Mb  ({loc['label']})")
    ax1.legend(loc="upper right", fontsize=8)
    ax1.set_xlim(lo/1e6, hi/1e6)

    # AAA panel
    ax2 = axes[1]
    ax2.scatter(aaa_sub["pos"]/1e6, aaa_sub["neg_log10P"],
                s=8, alpha=0.5, c="darkorange", label="AAA (AAAGen 2023)")
    aaa_lead_p = aaa_sub["P"].min()
    aaa_lead_pos = aaa_sub.loc[aaa_sub["P"].idxmin(), "pos"]
    ax2.scatter([aaa_lead_pos/1e6], [-np.log10(aaa_lead_p)],
                s=80, c="darkred", marker="*", zorder=5, edgecolors="white")
    ax2.axhline(-np.log10(5e-8), color="red", ls="--", lw=0.8)
    ax2.set_ylabel("-log10(P)  [AAA]")
    ax2.set_xlabel(f"Position on chromosome {c} (Mb)")
    ax2.legend(loc="upper right", fontsize=8)
    ax2.set_xlim(lo/1e6, hi/1e6)

    # Gene track at bottom
    gene_sub = genes[(genes["chr"]==c) & (genes["end"]>=lo) & (genes["start"]<=hi)]
    for _, g in gene_sub.iterrows():
        gx = max(g["start"], lo)/1e6
        gw = (min(g["end"], hi) - max(g["start"], lo))/1e6
        ax2.add_patch(Rectangle((gx, -2.5), gw, 1.5, color="gray", alpha=0.5))
        ax2.text(gx + gw/2, -1.7, g["gene"], ha="center", va="bottom",
                 fontsize=6, rotation=45)
    ax2.set_ylim(bottom=-3)

    plt.tight_layout()
    out = os.path.join(RESULTS_DIR, f"locus_plot_chr{c}_{loc['center']//1000}.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}  (IA pts={len(ia_sub)}, AAA pts={len(aaa_sub)})")
    print(f"    IA lead: pos={ia_lead_pos}, P={ia_lead_p:.2e}")
    print(f"    AAA lead: pos={aaa_lead_pos}, P={aaa_lead_p:.2e}")

for loc in TOP_LOCI:
    plot_locus(loc)

print("Done.")
