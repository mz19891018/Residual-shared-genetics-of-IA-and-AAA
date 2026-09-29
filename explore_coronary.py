import pyarrow.parquet as pq

fp = r'D:\IA_AAA_analysis\data\gtex\GTEx_Analysis_v11_eQTL\Artery_Coronary.v11.eQTLs.signif_pairs.parquet'
pf = pq.ParquetFile(fp)
df = pf.read().to_pandas()
print('Artery_Coronary:')
print('  pairs: %d' % len(df))
print('  genes: %d' % df['phenotype_id'].nunique())
print('  variants: %d' % df['variant_id'].nunique())
