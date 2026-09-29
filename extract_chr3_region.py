#!/usr/bin/env python
"""Extract chr3:6.4-7.5Mb region from original IA and AAA GWAS, same format as existing region files."""
import pandas as pd
import numpy as np
import os

BASE = r"D:\IA_AAA_analysis"
IA_FILE = os.path.join(BASE, "data", "gwas", "ia_2020", "IA.GWAS.BakkerMK.2020.sumstats.Stage_1.txt.gz")
AAA_FILE = os.path.join(BASE, "data", "gwas", "aaa_2023", "meta-AAAgen-final-sumstat.txt.gz")
OUT_DIR = os.path.join(BASE, "results", "susie")

CHROM = "3"
START = 6400000
END = 7500000
LOCUS = "chr3q12"

# ---- IA ----
print("Loading IA GWAS...")
ia = pd.read_csv(IA_FILE, sep=r'\s+', compression='gzip',
                 usecols=['CHR', 'BP', 'SNP', 'A_EFF', 'A_NONEFF', 'Freq_EFF', 'BETA', 'SE', 'P', 'MAF', 'Neff'])
ia = ia.rename(columns={'CHR':'chr','BP':'pos','SNP':'snp','A_EFF':'a1','A_NONEFF':'a2',
                        'Freq_EFF':'frq','BETA':'beta','SE':'se','P':'p','MAF':'maf','Neff':'n'})
ia['chr'] = ia['chr'].astype(str)
ia_reg = ia[(ia['chr']==CHROM) & (ia['pos']>=START) & (ia['pos']<=END)].copy()
ia_reg['locus'] = LOCUS
ia_out = os.path.join(OUT_DIR, f"region_chr3_{START}_{END}_IA.tsv")
ia_reg.to_csv(ia_out, sep='\t', index=False)
print(f"IA chr3 region: {len(ia_reg)} SNPs -> {ia_out}")
del ia, ia_reg

# ---- AAA ----
print("Loading AAA GWAS...")
aaa = pd.read_csv(AAA_FILE, sep=r'\s+', compression='gzip',
                  usecols=['chr','pos','ref','alt','N','AF','Effectsize','Effectsize_SD','pvalue'])
aaa = aaa.rename(columns={'chr':'chr','pos':'pos','ref':'a2','alt':'a1','AF':'frq',
                          'Effectsize':'beta','Effectsize_SD':'se','pvalue':'p','N':'n'})
aaa['chr'] = aaa['chr'].astype(str)
aaa['snp'] = aaa['chr'] + ":" + aaa['pos'].astype(str) + ":" + aaa['a2'] + ":" + aaa['a1']
aaa['maf'] = np.minimum(aaa['frq'], 1-aaa['frq'])
aaa_reg = aaa[(aaa['chr']==CHROM) & (aaa['pos']>=START) & (aaa['pos']<=END)].copy()
aaa_reg['locus'] = LOCUS
aaa_out = os.path.join(OUT_DIR, f"region_chr3_{START}_{END}_AAA.tsv")
aaa_reg.to_csv(aaa_out, sep='\t', index=False)
print(f"AAA chr3 region: {len(aaa_reg)} SNPs -> {aaa_out}")
del aaa, aaa_reg

print("Done.")
