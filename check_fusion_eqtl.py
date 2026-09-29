import pandas as pd

# Check eQTL counts
print('=== eQTL counts ===')
for f in ['MTAP_Artery_Tibial_v8EUR_eQTL', 'MTAP_Artery_Aorta_v8EUR_eQTL', 'CDKN2A_Artery_Aorta_v8EUR_eQTL']:
    try:
        df = pd.read_csv(rf'D:\IA_AAA_analysis\results\coloc\{f}.tsv', sep='\t')
        print(f'{f}: {len(df)} eQTLs, min p={df["pval_nominal"].min():.2e}')
    except Exception as e:
        print(f'{f}: ERROR {e}')

# Check FUSION results for 9p21.3
print('\n=== FUSION IA chr9 9p21.3 ===')
df_ia = pd.read_csv(r'D:\IA_AAA_analysis\results\twas\fusion\IA_chr9_aorta', sep='\t')
p21_ia = df_ia[(df_ia['CHR']==9) & (df_ia['P0']>=21000000) & (df_ia['P0']<=23000000)]
print(p21_ia[['ID','P0','TWAS.Z','TWAS.P']].to_string())

print('\n=== FUSION AAA chr9 9p21.3 ===')
df_aaa = pd.read_csv(r'D:\IA_AAA_analysis\results\twas\fusion\AAA_chr9_aorta', sep='\t')
p21_aaa = df_aaa[(df_aaa['CHR']==9) & (df_aaa['P0']>=21000000) & (df_aaa['P0']<=23000000)]
print(p21_aaa[['ID','P0','TWAS.Z','TWAS.P']].to_string())

# Check all significant genes in chr9
print('\n=== FUSION IA chr9 significant (P<0.05) ===')
sig_ia = df_ia[df_ia['TWAS.P'] < 0.05].sort_values('TWAS.P')
print(f'Total significant: {len(sig_ia)}')
print(sig_ia[['ID','P0','TWAS.Z','TWAS.P']].head(10).to_string())

print('\n=== FUSION AAA chr9 significant (P<0.05) ===')
sig_aaa = df_aaa[df_aaa['TWAS.P'] < 0.05].sort_values('TWAS.P')
print(f'Total significant: {len(sig_aaa)}')
print(sig_aaa[['ID','P0','TWAS.Z','TWAS.P']].head(10).to_string())
