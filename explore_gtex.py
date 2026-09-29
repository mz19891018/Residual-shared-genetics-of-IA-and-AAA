import pyarrow.parquet as pq

tissues = ['Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial']
results = {}
for t in tissues:
    fp = r'D:\IA_AAA_analysis\data\gtex\GTEx_Analysis_v11_eQTL\%s.v11.eQTLs.signif_pairs.parquet' % t
    pf = pq.ParquetFile(fp)
    n_rows = pf.metadata.num_rows
    # Read all data for unique counts
    df = pf.read().to_pandas()
    n_genes = df['phenotype_id'].nunique()
    n_vars = df['variant_id'].nunique()
    results[t] = (n_rows, n_genes, n_vars)
    print('%s:' % t)
    print('  Significant variant-gene pairs: %s' % format(n_rows, ','))
    print('  Unique genes (phenotype_id): %s' % format(n_genes, ','))
    print('  Unique variants: %s' % format(n_vars, ','))
    print('  Columns: %s' % list(df.columns))
    print()
