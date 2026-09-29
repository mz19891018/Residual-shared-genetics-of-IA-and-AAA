#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Comprehensive revision of Manuscript v3 to fix 12 issues identified by reviewer.
"""
from docx import Document
from docx.shared import Pt
from docx.oxml.ns import qn
import copy

doc = Document(r'D:\IA_AAA_analysis\Manuscript_v3_final.docx')

def set_paragraph_text(para, new_text):
    """Replace entire paragraph text while keeping first run's formatting"""
    if para.runs:
        para.runs[0].text = new_text
        for run in para.runs[1:]:
            run.text = ''
    else:
        para.add_run(new_text)

def remove_bold(para):
    """Remove bold from all runs"""
    for run in para.runs:
        run.bold = False

# ============================================================
# 1. Fix Abstract Results paragraph (Para 9) - remove bold, deduplicate, shorten
# ============================================================
abstract_results = (
    "Results: IA and AAA were genetically correlated (rg=0.35, SE 0.06, P=5.0×10⁻⁸) "
    "without evidence of sample overlap. Smoking initiation was causally associated with both diseases. "
    "SBP was causally associated with IA; for AAA, conventional IVW was null (P=0.20) but outlier-robust methods "
    "(MR-RAPS, MR-PRESSO) supported a moderate causal effect (β≈0.010–0.012 per mmHg, P<10⁻¹⁰), "
    "though heterogeneity was substantial (Q=1676). After joint conditioning on SBP and smoking, "
    "the residual correlation was 0.26 (95% CI 0.10–0.42, P=0.0014) and remained essentially unchanged "
    "after further conditioning on LDL cholesterol (0.27) or CAD (0.25). Local correlation was strongest "
    "at 9p21.3 (ρ=0.73, P=4.0×10⁻¹³) and 18q11 (ρ=0.68, P=4.9×10⁻⁷), but no region remained significant "
    "after conditioning. Colocalization supported shared causal variants at 9p21.3 (PP.H4=0.88) and "
    "18q11 (PP.H4=0.998). Neither signal colocalized with an arterial eQTL: CDKN2A had few arterial eQTLs, "
    "and MTAP, the dominant arterial eGene at 9p21.3, did not colocalize with either disease."
)

for i, para in enumerate(doc.paragraphs):
    if para.text.startswith('SNP heritability was significant') or (
        'Results: IA and AAA were genetically correlated' in para.text and i < 15):
        remove_bold(para)
        set_paragraph_text(para, abstract_results)
        print(f'Fixed Para {i}: Abstract Results rewritten (~{len(abstract_results.split())} words)')
        break

# ============================================================
# 2. Restore eQTL section body (Para 52 - currently empty)
# ============================================================
eqtl_text = (
    "We tested colocalization of the shared 9p21.3 and 18q11 signals with cis-eQTLs for all genes "
    "within ±1 Mb in GTEx artery tissues (tibial, aorta, coronary). At 9p21.3, CDKN2A had only one "
    "significant eQTL variant in tibial artery and one in aorta, and was not a strong eGene in arterial "
    "tissue. MTAP was the principal eGene at 9p21.3 (530 significant eQTL variants in tibial artery; "
    "minimum P=1.1×10⁻²¹), but its eQTL signal did not colocalize with IA or AAA (z-correlation r=0.06–0.13, "
    "all P>0.05); the MTAP eQTLs lie about 200 kb proximal to the aneurysm association peak. CDKN2B had "
    "no significant arterial eQTLs. An earlier analysis suggesting colocalization with CDKN2A rested on "
    "very few variants and was not robust. In LD-aware TWAS of chromosome 9, no 9p21.3 gene was associated "
    "with IA (MTAP P=0.62) or with AAA after correction for the number of genes tested (MTAP P=0.063). "
    "No gene colocalized with the 18q11 signal. Formal eQTL colocalization with full summary statistics "
    "from eQTL Catalogue is planned as future work."
)

for i, para in enumerate(doc.paragraphs):
    if para.style.name == 'Heading 3' and 'No arterial eQTL' in para.text:
        # Next paragraph should be the body
        if i+1 < len(doc.paragraphs) and not doc.paragraphs[i+1].text.strip():
            set_paragraph_text(doc.paragraphs[i+1], eqtl_text)
            print(f'Fixed Para {i+1}: Restored eQTL section body')
        break

# ============================================================
# 3. Fix TWAS paragraph (Para 58) - remove fragmented sentences
# ============================================================
twas_text = (
    "In an exploratory TWAS using FUSION arterial expression weights [37], a rapid scan identified "
    "20 genes nominally associated with both IA and AAA (Supplementary Table S2), mapping to about 15 "
    "independent loci; for example, GBA2, NPR2 and RGP1 all lie at 9p13.3. Fourteen showed concordant "
    "directions of effect in the two diseases and six showed opposite directions. However, this scan did "
    "not account for LD and likely inflated Z-statistics. FUSION LD-corrected analysis of chromosome 9 "
    "(248 genes tested) showed no significant genes at 9p21.3 after Bonferroni correction (threshold "
    "P≈2×10⁻⁴); the smallest P-value was 0.007 for an unannotated lncRNA (ENSG00000264801.1) in AAA. "
    "CDKN2A did not pass FUSION's model-performance filter and was therefore not tested. Full-genome "
    "FUSION LD-corrected analysis is ongoing."
)

for i, para in enumerate(doc.paragraphs):
    if 'exploratory TWAS using FUSION' in para.text:
        set_paragraph_text(para, twas_text)
        print(f'Fixed Para {i}: TWAS paragraph rewritten')
        break

# ============================================================
# 4. Fix single-cell paragraph (Para 54) - unify comparison description
# ============================================================
sc_text = (
    "In the murine IA model (GSE193533; Formed vs Sham aneurysms, n=1 pooled sample per condition, "
    "no biological replicates), Cdkn2a expression was markedly higher in synthetic-like than in contractile "
    "VSMCs (about 36-fold) and was highest in fibroblasts (Figure 5). Because each condition consisted of "
    "a single pooled sample, these differences are descriptive. In human AAA tissue (GSE166676), CDKN2A "
    "was detected in a higher proportion of synthetic-like than contractile VSMCs, but VSMCs were few and "
    "cells from AAA and control samples clustered separately, which precluded formal comparison. RBBP8 was "
    "expressed at low levels without state specificity, and GATA6 and CHRNA5 were barely detected in "
    "vessel-wall cells."
)

for i, para in enumerate(doc.paragraphs):
    if 'murine IA model, Cdkn2a expression was markedly higher' in para.text:
        set_paragraph_text(para, sc_text)
        print(f'Fixed Para {i}: Single-cell paragraph rewritten')
        break

# ============================================================
# 5. Fix Discussion first paragraph (Para 60) - CDKN2A as candidate not conclusion
# ============================================================
for i, para in enumerate(doc.paragraphs):
    if 'shared 9p21.3 signal pointed to CDKN2A' in para.text:
        new_text = para.text.replace(
            'the shared 9p21.3 signal pointed to CDKN2A and smooth muscle cell state as candidate mechanisms.',
            'the shared 9p21.3 signal did not colocalize with an arterial eQTL; CDKN2A and smooth muscle cell state remain candidate mechanisms requiring functional validation.'
        )
        set_paragraph_text(para, new_text)
        print(f'Fixed Para {i}: Discussion CDKN2A wording')
        break

# ============================================================
# 6. Fix Supplementary Table S1 footnote - clarify Egger intercept
# ============================================================
for i, para in enumerate(doc.paragraphs):
    if 'MR-Egger intercept was not significant (P=0.76)' in para.text:
        new_text = para.text.replace(
            'the MR-Egger intercept was not significant (P=0.76), arguing against directional pleiotropy.',
            'the MR-Egger intercept was not significant in the robust analysis (P=0.76), arguing against directional pleiotropy; the conventional MR-Egger intercept P for SBP→AAA was 0.080 (Table S1).'
        )
        set_paragraph_text(para, new_text)
        print(f'Fixed Para {i}: Egger intercept clarification')
        break

# ============================================================
# 7. Fix Table 4 22q12 candidate genes - remove distant genes
# ============================================================
for table in doc.tables:
    for row in table.rows:
        if row.cells[0].text.strip() == '22q12':
            for cell in row.cells:
                if 'LIMK2' in cell.text:
                    for para in cell.paragraphs:
                        set_paragraph_text(para, 'LIMK2')
                    print('Fixed Table 4: 22q12 candidate genes -> LIMK2 only')

# ============================================================
# 8. Fix Table 1 height reference - Yengo 2018 not 2022
# ============================================================
for table in doc.tables:
    for row in table.rows:
        if 'Height' in row.cells[0].text:
            for cell in row.cells:
                if 'Yengo et al. 2022' in cell.text:
                    for para in cell.paragraphs:
                        set_paragraph_text(para, 'GIANT, Yengo et al. 2018 [14]')
                    print('Fixed Table 1: Height reference -> Yengo 2018')

# ============================================================
# 9. Fix references [31] and [32] - GEO original papers
# ============================================================
ref_31 = "Davis FM, Tsoi LC, Melvin WJ, denDekker A, et al. Inhibition of macrophage histone demethylase JMJD3 protects against abdominal aortic aneurysms. J Exp Med. 2021;218(6). doi:10.1084/jem.20201323."
ref_32 = "Martinez G, et al. Single-cell transcriptome analysis of the circle of Willis in a mouse cerebral aneurysm model. Stroke. 2022;53:2647–2657. doi:10.1161/STROKEAHA.121.037126."

ref_count = 0
for i, para in enumerate(doc.paragraphs):
    if 'Single-cell analysis of abdominal aortic aneurysm reveals vascular smooth muscle cell heterogeneity' in para.text:
        set_paragraph_text(para, ref_31)
        print(f'Fixed Para {i}: Reference [31] -> Davis FM 2021 J Exp Med')
        ref_count += 1
    elif 'Single-cell transcriptomics of intracranial aneurysm in a murine model' in para.text:
        set_paragraph_text(para, ref_32)
        print(f'Fixed Para {i}: Reference [32] -> Martinez 2022 Stroke')
        ref_count += 1

# ============================================================
# 10. Add numbering to references
# ============================================================
ref_start = None
for i, para in enumerate(doc.paragraphs):
    if para.style.name == 'Heading 2' and 'References' in para.text:
        ref_start = i + 1
        break

if ref_start:
    ref_num = 1
    for i in range(ref_start, len(doc.paragraphs)):
        para = doc.paragraphs[i]
        text = para.text.strip()
        if text and not text.startswith('['):
            # Add number
            new_text = f'[{ref_num}] {text}'
            set_paragraph_text(para, new_text)
            ref_num += 1
    print(f'Added numbering to {ref_num - 1} references')

# ============================================================
# 11. Fix height reference in reference list [14]
# ============================================================
for i, para in enumerate(doc.paragraphs):
    if 'Yengo L, Vedantam S, Marouli E' in para.text and '2022' in para.text:
        new_ref = "[14] Yengo L, Sidorenko J, Kemper KE, et al. Meta-analysis of genome-wide association studies for height and body mass index in ~700000 individuals of European ancestry. Hum Mol Genet. 2018;27(20):3641–3649. doi:10.1093/hmg/ddy271."
        set_paragraph_text(para, new_ref)
        print(f'Fixed Para {i}: Reference [14] -> Yengo 2018 Hum Mol Genet')
        break

# ============================================================
# 12. Check for orphan punctuation/spaces
# ============================================================
for i, para in enumerate(doc.paragraphs):
    text = para.text
    # Fix double spaces
    if '  ' in text:
        for run in para.runs:
            if '  ' in run.text:
                run.text = run.text.replace('  ', ' ')
    # Fix " ; " patterns
    if ' ; ' in text:
        for run in para.runs:
            if ' ; ' in run.text:
                run.text = run.text.replace(' ; ', '; ')

print('\n=== All fixes applied ===')

# Save
output_path = r'D:\IA_AAA_analysis\Manuscript_v4.docx'
doc.save(output_path)
print(f'Saved to: {output_path}')
print(f'Paragraphs: {len(doc.paragraphs)}, Tables: {len(doc.tables)}')
