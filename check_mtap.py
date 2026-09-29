import pandas as pd

for tissue in ['Tibial', 'Aorta']:
    df = pd.read_csv(rf'D:\IA_AAA_analysis\results\coloc\MTAP_Artery_{tissue}_v11_eQTL.tsv', sep='\t')
    print(f'MTAP {tissue}: {len(df)} eQTLs')
    print(f'  pval range: {df["pval_nominal"].min():.2e} to {df["pval_nominal"].max():.2e}')
    print(f'  slope range: {df["slope"].min():.3f} to {df["slope"].max():.3f}')
    print(f'  variant_id example: {df["variant_id"].iloc[0]}')

# Check 9p21.3 GWAS
for trait in ['IA_2020', 'AAA_2023']:
    df = pd.read_csv(rf'D:\IA_AAA_analysis\results\coloc\{trait}_9p21_3.tsv', sep='\t')
    print(f'{trait} 9p21.3: {len(df)} SNPs, min P={df["P"].min():.2e}')
