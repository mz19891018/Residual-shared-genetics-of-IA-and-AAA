import gzip
import pandas as pd

# Correct IDs: MTAP=ENSG00000099810, CDKN2A=ENSG00000147889
genes = {
    'ENSG00000099810': 'MTAP',
    'ENSG00000147889': 'CDKN2A',
    'ENSG00000147873': 'CDKN2B',
}

for tissue in ['Artery_Tibial', 'Artery_Aorta']:
    print(f'\n=== {tissue} ===')
    inp = rf'D:\IA_AAA_analysis\data\gtex\GTEx_Analysis_v8_eQTL_EUR\eqtls\{tissue}.v8.EUR.signif_pairs.txt.gz'
    with gzip.open(inp, 'rt') as f:
        header = f.readline().strip().split('\t')
        pid_idx = header.index('phenotype_id')
        counts = {g: 0 for g in genes}
        rows = {g: [] for g in genes}
        for line in f:
            parts = line.strip().split('\t')
            pid = parts[pid_idx]
            for ensg in genes:
                if pid.startswith(ensg):
                    counts[ensg] += 1
                    rows[ensg].append(parts)
                    break
        for ensg, name in genes.items():
            print(f'  {name} ({ensg}): {counts[ensg]} significant eQTL pairs')
            if counts[ensg] > 0:
                # Save
                out = rf'D:\IA_AAA_analysis\results\coloc\{name}_{tissue}_v8EUR_eQTL.tsv'
                with open(out, 'w') as o:
                    o.write('\t'.join(header) + '\n')
                    for r in rows[ensg]:
                        o.write('\t'.join(r) + '\n')
                print(f'    Saved to {out}')
                # Show top 3
                df = pd.DataFrame(rows[ensg], columns=header)
                df['pval_nominal'] = df['pval_nominal'].astype(float)
                df = df.sort_values('pval_nominal')
                print(df[['variant_id','pval_nominal','slope','slope_se']].head(3).to_string())
