import pyarrow.parquet as pq
import pyarrow as pa

# Check v11 Artery Tibial for key genes in 9p21.3
genes = {
    'ENSG00000147889': 'CDKN2A',
    'ENSG00000147883': 'CDKN2B',  # verify
    'ENSG00000240498': 'CDKN2B-AS1',  # ANRIL, verify
    'ENSG00000134184': 'MTAP',  # verify
}

# Read with filter for efficiency
print('=== GTEx v11 Artery Tibial ===')
pf = pq.ParquetFile(r'D:\IA_AAA_analysis\data\gtex\GTEx_Analysis_v11_eQTL\Artery_Tibial.v11.eQTLs.signif_pairs.parquet')
print(f'Row groups: {pf.metadata.num_row_groups}, rows: {pf.metadata.num_rows}')

# Read in chunks and filter
results = {g: [] for g in genes}
for batch in pf.iter_batches(batch_size=100000):
    df = batch.to_pandas()
    for ensg, name in genes.items():
        mask = df['phenotype_id'].astype(str).str.contains(ensg, na=False)
        if mask.any():
            results[ensg].extend(df[mask].to_dict('records'))

for ensg, name in genes.items():
    print(f'{name} ({ensg}): {len(results[ensg])} eQTL pairs')
    if len(results[ensg]) > 0:
        r = results[ensg][0]
        print(f'  Example: variant={r.get("variant_id")}, pval={r.get("pval_nominal")}, slope={r.get("slope")}')

print('\n=== GTEx v11 Artery Aorta ===')
pf2 = pq.ParquetFile(r'D:\IA_AAA_analysis\data\gtex\GTEx_Analysis_v11_eQTL\Artery_Aorta.v11.eQTLs.signif_pairs.parquet')
results2 = {g: [] for g in genes}
for batch in pf2.iter_batches(batch_size=100000):
    df = batch.to_pandas()
    for ensg, name in genes.items():
        mask = df['phenotype_id'].astype(str).str.contains(ensg, na=False)
        if mask.any():
            results2[ensg].extend(df[mask].to_dict('records'))

for ensg, name in genes.items():
    print(f'{name} ({ensg}): {len(results2[ensg])} eQTL pairs')
