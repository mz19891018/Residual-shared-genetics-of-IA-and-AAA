#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Final fixes: add heritability results, fix empty parentheses, remove reference note"""
from docx import Document

doc = Document(r'D:\IA_AAA_analysis\Manuscript_v3_final.docx')

# 1. Add heritability results after Para 38 (genetic correlation result)
# Find the paragraph and insert heritability info
for i, para in enumerate(doc.paragraphs):
    if 'IA and AAA were genetically correlated' in para.text:
        # Add heritability to this paragraph
        full_text = para.text
        if 'heritab' not in full_text.lower():
            # Insert at beginning
            new_text = ('SNP heritability was significant for both traits '
                       '(IA h²=0.20, SE=0.03; AAA h²=0.014, SE=0.001 on the observed scale; '
                       'liability-scale estimates assuming K=0.03 for IA and K=0.02 for AAA). '
                       + full_text)
            # Set text in first run, clear others
            if para.runs:
                para.runs[0].text = new_text
                for run in para.runs[1:]:
                    run.text = ''
            print(f'Fixed Para {i}: added heritability results')
        break

# 2. Fix empty parentheses in Para 39
for i, para in enumerate(doc.paragraphs):
    if 'AAA was not (rg=0.05, SE=0.03, P=0.053; )' in para.text:
        for run in para.runs:
            if 'P=0.053; )' in run.text:
                run.text = run.text.replace('P=0.053; )', 'P=0.053)')
                print(f'Fixed Para {i}: removed empty parentheses')
        break

# 3. Remove reference note (Para 90)
for i, para in enumerate(doc.paragraphs):
    if '卷期页码请用文献管理软件核对' in para.text:
        # Clear the paragraph
        for run in para.runs:
            run.text = ''
        print(f'Fixed Para {i}: removed reference note')
        break

# 4. Check for any remaining issues
print('\n=== Final check ===')
full_text = '\n'.join([p.text for p in doc.paragraphs])

# Check heritability
if 'h²=0.20' in full_text or 'h2=0.20' in full_text:
    print('  [OK] IA heritability reported')
else:
    print('  [MISSING] IA heritability')

if '0.014' in full_text:
    print('  [OK] AAA heritability reported')
else:
    print('  [MISSING] AAA heritability')

# Check no empty parentheses
if '; )' in full_text:
    print('  [WARNING] Empty parentheses remain')
else:
    print('  [OK] No empty parentheses')

# Check no 【待核对】 markers
if '【待核对' in full_text or '【待重算' in full_text:
    print('  [WARNING] 待核对 markers remain')
else:
    print('  [OK] No 待核对 markers')

# Save
output_path = r'D:\IA_AAA_analysis\Manuscript_v3_final.docx'
doc.save(output_path)
print(f'\nSaved to: {output_path}')
print(f'Paragraphs: {len(doc.paragraphs)}, Tables: {len(doc.tables)}')
