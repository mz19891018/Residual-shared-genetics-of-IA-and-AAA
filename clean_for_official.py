#!/usr/bin/env python3
"""Clean munged sumstats for official LDSC:
- keep only biallelic single-base SNPs (A1,A2 in {A,C,G,T})
- keep only SNPs present in HapMap3 no-MHC LD reference list
This is equivalent to official munge_sumstats.py --merge-alleles.
"""
import pandas as pd, os, gzip

QC = r'D:\IA_AAA_analysis\data\qc'
OUT = r'D:\IA_AAA_analysis\data\qc\official'
os.makedirs(OUT, exist_ok=True)

# HapMap3 no-MHC SNP list used by LD reference
snplist = pd.read_csv(r'D:\IA_AAA_analysis\data\ld_reference\hm3_no_MHC.list.txt',
                      sep=r'\s+', header=None, names=['SNP'])
hm3 = set(snplist.SNP)
print(f'HapMap3 SNP list: {len(hm3)} SNPs')

TRAITS = [
    'IA_2020', 'AAA_2023', 'SBP', 'SmokingInit', 'CAD_2022',
    'LDL', 'HDL', 'TG', 'Height',
    'AA_diam_P2020', 'AscAortaDiam_P2023', 'AscAortaDistens_P2023',
    'AscAortaStrain_P2023', 'DA_diam_P2020', 'DescAortaDiam_P2023',
    'DescAortaDistens_P2023', 'DescAortaStrain_P2023',
]

BASES = {'A', 'C', 'G', 'T'}
for t in TRAITS:
    src = os.path.join(QC, t + '.sumstats.gz')
    if not os.path.exists(src):
        print(f'MISSING: {src}')
        continue
    df = pd.read_csv(src, sep='\t')
    n0 = len(df)
    # biallelic single-base only
    df = df[df.A1.isin(BASES) & df.A2.isin(BASES)]
    # restrict to HapMap3
    df = df[df.SNP.isin(hm3)]
    # drop duplicates on SNP
    df = df.drop_duplicates(subset='SNP', keep='first')
    dst = os.path.join(OUT, t + '.sumstats.gz')
    df.to_csv(dst, sep='\t', index=False, compression='gzip')
    print(f'{t}: {n0} -> {len(df)} SNPs (dropped {n0-len(df)})')
print('DONE')
