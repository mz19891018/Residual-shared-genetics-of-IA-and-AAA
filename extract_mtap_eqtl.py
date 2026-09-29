import pyarrow.parquet as pq
import pandas as pd
import gzip

# Extract MTAP eQTLs from GTEx v11 Artery Tibial
print('Extracting MTAP eQTLs from GTEx v11 Artery Tibial...')
pf = pq.ParquetFile(r'D:\IA_AAA_analysis\data\gtex\GTEx_Analysis_v11_eQTL\Artery_Tibial.v11.eQTLs.signif_pairs.parquet')
mtap_eqtl = []
for batch in pf.iter_batches(batch_size=100000):
    df = batch.to_pandas()
    mask = df['phenotype_id'].astype(str).str.contains('ENSG00000134184', na=False)
    if mask.any():
        mtap_eqtl.append(df[mask])

if mtap_eqtl:
    mtap_df = pd.concat(mtap_eqtl)
    print(f'MTAP eQTLs: {len(mtap_df)}')
    print(f'Columns: {list(mtap_df.columns)}')
    # Save for colocalization
    mtap_df.to_csv(r'D:\IA_AAA_analysis\results\coloc\MTAP_Artery_Tibial_v11_eQTL.tsv', sep='\t', index=False)
    print('Saved MTAP eQTLs')
    print(mtap_df[['variant_id','pval_nominal','slope','slope_se']].head(5).to_string())
else:
    print('No MTAP eQTLs found')

# Also extract from Artery Aorta
print('\nExtracting MTAP eQTLs from GTEx v11 Artery Aorta...')
pf2 = pq.ParquetFile(r'D:\IA_AAA_analysis\data\gtex\GTEx_Analysis_v11_eQTL\Artery_Aorta.v11.eQTLs.signif_pairs.parquet')
mtap_eqtl2 = []
for batch in pf2.iter_batches(batch_size=100000):
    df = batch.to_pandas()
    mask = df['phenotype_id'].astype(str).str.contains('ENSG00000134184', na=False)
    if mask.any():
        mtap_eqtl2.append(df[mask])

if mtap_eqtl2:
    mtap_df2 = pd.concat(mtap_eqtl2)
    print(f'MTAP eQTLs in Aorta: {len(mtap_df2)}')
    mtap_df2.to_csv(r'D:\IA_AAA_analysis\results\coloc\MTAP_Artery_Aorta_v11_eQTL.tsv', sep='\t', index=False)
else:
    print('No MTAP eQTLs in Aorta')

# Extract 9p21.3 GWAS data for IA and AAA
print('\nExtracting 9p21.3 GWAS data...')
for trait in ['IA_2020', 'AAA_2023']:
    inp = rf'D:\IA_AAA_analysis\data\qc\{trait}.sumstats.gz'
    out = rf'D:\IA_AAA_analysis\results\coloc\{trait}_9p21_3.tsv'
    # Need chr/bp - munged file doesn't have it, use 1000G bim to map
    # First read 1000G bim for chr9
    bim = pd.read_csv(r'D:\IA_AAA_analysis\data\ld_reference\1000G_EUR_Phase3_plink\1000G.EUR.QC.9.bim',
                      sep='\t', header=None, names=['CHR','SNP','CM','BP','A1','A2'])
    chr9_snps = set(bim[(bim['BP']>=21500000) & (bim['BP']<=22500000)]['SNP'])
    
    with gzip.open(inp, 'rt') as f, open(out, 'w') as o:
        header = f.readline()
        o.write(header)
        count = 0
        for line in f:
            snp = line.split('\t')[0]
            if snp in chr9_snps:
                o.write(line)
                count += 1
    print(f'{trait}: {count} SNPs in 9p21.3')
