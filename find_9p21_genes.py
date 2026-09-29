import gzip
import pandas as pd

# Find genes in 9p21.3 region from GTEx v8 egenes
print('=== Genes in 9p21.3 (chr9:21-23Mb) from GTEx v8 Artery Tibial egenes ===')
with gzip.open(r'D:\IA_AAA_analysis\data\gtex\GTEx_Analysis_v8_eQTL_EUR\eqtls\Artery_Tibial.v8.EUR.egenes.txt.gz', 'rt') as f:
    header = f.readline().strip().split('\t')
    print(f'Columns: {header[:5]}...')
    vid_idx = header.index('variant_id')
    pid_idx = header.index('phenotype_id')
    pval_idx = header.index('pval_nominal')
    for line in f:
        parts = line.strip().split('\t')
        vid = parts[vid_idx]
        # variant_id format: chr9_21999xxx_A_G_b38
        if vid.startswith('chr9_'):
            pos = int(vid.split('_')[1])
            if 21000000 <= pos <= 23000000:
                print(f'  {parts[pid_idx]}  lead_variant={vid}  pval={parts[pval_idx]}')

print('\n=== Searching signif_pairs for chr9:21-22Mb eQTLs ===')
count = 0
genes_in_region = set()
with gzip.open(r'D:\IA_AAA_analysis\data\gtex\GTEx_Analysis_v8_eQTL_EUR\eqtls\Artery_Tibial.v8.EUR.signif_pairs.txt.gz', 'rt') as f:
    header = f.readline().strip().split('\t')
    print(f'signif_pairs columns: {header}')
    vid_idx = header.index('variant_id')
    pid_idx = header.index('phenotype_id')
    for line in f:
        parts = line.strip().split('\t')
        vid = parts[vid_idx]
        if vid.startswith('chr9_'):
            pos = int(vid.split('_')[1])
            if 21500000 <= pos <= 22500000:
                genes_in_region.add(parts[pid_idx])
                count += 1
                if count <= 10:
                    print(f'  {parts[pid_idx]}  {vid}  pval={parts[header.index("pval_nominal")]}')

print(f'\nTotal eQTL pairs in 9p21.3: {count}')
print(f'Unique genes in 9p21.3: {len(genes_in_region)}')
for g in sorted(genes_in_region):
    print(f'  {g}')
