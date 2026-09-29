"""
CAD overlap check: For each of the 7 IA-AAA overlapping loci,
find the minimum CAD P-value within the region.
Uses LD score files to map rsID -> (chr, bp), then looks up CAD P.
"""
import gzip
import csv
import os

# Paths
LD_DIR = r"D:\IA_AAA_analysis\data\ld_reference\LDscore"
CAD_FILE = r"D:\IA_AAA_analysis\data\qc\CAD_2022.sumstats.gz"
OUT_FILE = r"D:\IA_AAA_analysis\results\cad_overlap_check.tsv"

# 7 overlapping loci from locus_overlap.tsv
# (locus_id, chr, ov_start, ov_end, ia_lead_snp, ia_lead_pos, ia_lead_p)
LOCI = [
    (1, 9,  21603341, 22598574, "rs1537373",   22103341, 2.595e-22),
    (2, 18, 19782856, 20723695, "rs11661542",  20223695, 5.735e-16),
    (3, 15, 78367482, 79314046, "rs10519203",  78814046, 1.293e-09),
    (4, 16, 74951461, 75812548, "rs11646044",  75312548, 5.209e-11),
    (5, 10, 105140978,105452499,"rs79780963",  104952499, 6.817e-09),
    (6, 22, 30498650, 30843186, "rs39713",     30343186, 4.101e-08),
    (7, 22, 29843186, 30498651, "rs39713",     30343186, 4.101e-08),
]

print("Loading CAD munged sumstats...")
cad_p = {}  # rsID -> P
cad_z = {}
cad_beta = {}
cad_se = {}
with gzip.open(CAD_FILE, 'rt') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        snp = row['SNP']
        cad_p[snp] = float(row['P'])
        cad_z[snp] = float(row['Z'])
        cad_beta[snp] = float(row['BETA'])
        cad_se[snp] = float(row['SE'])
print(f"  Loaded {len(cad_p)} CAD SNPs")

print("Checking loci...")
results = []
for loc_id, chrom, start, end, lead_snp, lead_pos, lead_p in LOCI:
    # Load LD score for this chromosome to get SNP positions
    ld_file = os.path.join(LD_DIR, f"LDscore.{chrom}.l2.ldscore.gz")
    region_snps = []  # (bp, rsid)
    with gzip.open(ld_file, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            bp = int(row['BP'])
            if start <= bp <= end:
                region_snps.append((bp, row['SNP']))
    
    # Find min CAD P among region SNPs
    min_p = 1.0
    min_snp = None
    min_bp = None
    min_z = None
    min_beta = None
    min_se = None
    for bp, rsid in region_snps:
        if rsid in cad_p:
            p = cad_p[rsid]
            if p < min_p:
                min_p = p
                min_snp = rsid
                min_bp = bp
                min_z = cad_z[rsid]
                min_beta = cad_beta[rsid]
                min_se = cad_se[rsid]
    
    # Also check lead SNP directly
    lead_cad_p = cad_p.get(lead_snp, None)
    
    is_sig = "YES" if min_p < 5e-8 else "NO"
    classification = "atherosclerosis-related candidate" if min_p < 5e-8 else "non-atherosclerotic / IA-AAA specific"
    
    print(f"  Locus {loc_id} chr{chrom}:{start}-{end}: "
          f"region SNPs in CAD={sum(1 for _,r in region_snps if r in cad_p)}, "
          f"min P={min_p:.2e} ({min_snp}@{min_bp}), "
          f"lead {lead_snp} CAD P={lead_cad_p}, sig={is_sig}")
    
    results.append({
        'locus_id': loc_id,
        'chr': chrom,
        'region_start': start,
        'region_end': end,
        'ia_lead_snp': lead_snp,
        'ia_lead_pos': lead_pos,
        'ia_lead_p': lead_p,
        'cad_lead_snp_in_region': min_snp if min_snp else 'NA',
        'cad_lead_pos': min_bp if min_bp else 'NA',
        'cad_min_p': f"{min_p:.6e}",
        'cad_min_z': f"{min_z:.4f}" if min_z else 'NA',
        'cad_min_beta': f"{min_beta:.6f}" if min_beta else 'NA',
        'cad_min_se': f"{min_se:.6f}" if min_se else 'NA',
        'cad_lead_snp_p': f"{lead_cad_p:.6e}" if lead_cad_p else 'not_in_cad',
        'n_region_snps_in_cad': sum(1 for _, r in region_snps if r in cad_p),
        'cad_genomewide_sig': is_sig,
        'classification': classification,
    })

# Write output
fields = ['locus_id','chr','region_start','region_end','ia_lead_snp','ia_lead_pos','ia_lead_p',
          'cad_lead_snp_in_region','cad_lead_pos','cad_min_p','cad_min_z','cad_min_beta','cad_min_se',
          'cad_lead_snp_p','n_region_snps_in_cad','cad_genomewide_sig','classification']
with open(OUT_FILE, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
    w.writeheader()
    for r in results:
        w.writerow(r)
print(f"\nResults saved to {OUT_FILE}")
