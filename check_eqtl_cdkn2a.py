import pandas as pd
import pyarrow.parquet as pq

# GTEx v11 Artery Tibial
print('=== GTEx v11 Artery Tibial ===')
df = pq.read_table(r'D:\IA_AAA_analysis\data\gtex\GTEx_Analysis_v11_eQTL\Artery_Tibial.v11.eQTLs.signif_pairs.parquet').to_pandas()
print(f'Columns: {list(df.columns)}')
print(f'Total rows: {len(df)}')
cdkn2a = df[df['phenotype_id'].astype(str).str.contains('ENSG00000147889', na=False)]
print(f'CDKN2A eQTL pairs: {len(cdkn2a)}')
if len(cdkn2a) > 0:
    print(cdkn2a[['phenotype_id','variant_id','pval_nominal','slope','slope_se']].head(5).to_string())
    # Check variant_id format
    print(f'Variant ID example: {cdkn2a["variant_id"].iloc[0]}')

print('\n=== GTEx v11 Artery Aorta ===')
df2 = pq.read_table(r'D:\IA_AAA_analysis\data\gtex\GTEx_Analysis_v11_eQTL\Artery_Aorta.v11.eQTLs.signif_pairs.parquet').to_pandas()
cdkn2a2 = df2[df2['phenotype_id'].astype(str).str.contains('ENSG00000147889', na=False)]
print(f'CDKN2A eQTL pairs: {len(cdkn2a2)}')

print('\n=== GTEx v8 EUR Artery Tibial ===')
df3 = pd.read_csv(r'D:\IA_AAA_analysis\data\gtex\GTEx_Analysis_v8_eQTL_EUR\eqtls\Artery_Tibial.v8.EUR.signif_pairs.txt.gz', sep='\t')
print(f'Columns: {list(df3.columns)}')
gene_col = [c for c in df3.columns if 'gene' in c.lower() or 'phenotype' in c.lower()][0]
cdkn2a3 = df3[df3[gene_col].astype(str).str.contains('ENSG00000147889', na=False)]
print(f'CDKN2A eQTL pairs: {len(cdkn2a3)}')
if len(cdkn2a3) > 0:
    print(cdkn2a3.head(3).to_string())

print('\n=== GTEx v8 EUR Artery Aorta ===')
df4 = pd.read_csv(r'D:\IA_AAA_analysis\data\gtex\GTEx_Analysis_v8_eQTL_EUR\eqtls\Artery_Aorta.v8.EUR.signif_pairs.txt.gz', sep='\t')
cdkn2a4 = df4[df4[gene_col].astype(str).str.contains('ENSG00000147889', na=False)]
print(f'CDKN2A eQTL pairs: {len(cdkn2a4)}')
