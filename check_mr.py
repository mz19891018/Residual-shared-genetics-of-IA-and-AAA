#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd

# Check mr_summary.tsv
df = pd.read_csv(r'D:\IA_AAA_analysis\results\mr\mr_summary.tsv', sep='\t')
print('=== mr_summary.tsv columns ===')
print(df.columns.tolist())
print()
print('=== SBP to AAA rows ===')
for col in df.columns:
    if df[col].dtype == object:
        mask = df[col].str.contains('SBP', na=False) & df[col].str.contains('AAA', na=False)
        if mask.any():
            print(df[mask].to_string())
            break
print()

# Check if there's an egger intercept column
print('=== All columns ===')
for c in df.columns:
    print(f'  {c}')
