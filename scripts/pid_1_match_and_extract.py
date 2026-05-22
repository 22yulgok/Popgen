import os
import csv
import re
import numpy as np
import collections
import pandas as pd

BASE_DIR = "/home/duri_bae/popgen_lab"
OUT_DIR = os.path.join(BASE_DIR, "output/adaptive_sweep/PID_comparison")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, "figures"), exist_ok=True)

REPORT_FILE = os.path.join(BASE_DIR, "output/adaptive_sweep/adaptive_sweep_report_260522.md")
KNOWN_CATALOG = os.path.join(BASE_DIR, "output/adaptive_sweep/known_selection_snp_catalog_260522.csv")
CANDIDATES = os.path.join(BASE_DIR, "output/adaptive_sweep/high_frequency_polymorphic_candidates_260522.csv")

FAM_FILE = os.path.join(BASE_DIR, "raw/Tibetan.fstat.1240K.260508.fam")
BIM_FILE = os.path.join(BASE_DIR, "raw/Tibetan.fstat.1240K.260508.bim")
BED_FILE = os.path.join(BASE_DIR, "raw/Tibetan.fstat.1240K.260508.bed")

NEW_PIDS = ['Haneyi', 'Xingyi_EN', 'BGD_EN', 'YR_LN', 'Upper_YR_LN', 'Shandong_EN', 'Xiaojingshan', 'Yangshao_UYR', 'Shimao_group', 'STM_EN', 'Boisman_MN', 'Yumin', 'HMMH_MN', 'DevilsCave_N', 'Qihe3', 'Liangdao', 'Xitoucun', 'Xueshan']

# 1. Infer previous targets
targets = set()
with open(REPORT_FILE, 'r') as f:
    for line in f:
        if line.startswith("- **Target Populations Found**:"):
            pops_str = line.split(":", 1)[1].strip()
            for p in pops_str.split(","):
                p = p.strip()
                if p:
                    targets.add(p)
print(f"Inferred {len(targets)} previous targets from report.")

# 2. Match New PIDs in FAM
samples = []
pid_to_indices = collections.defaultdict(list)

with open(FAM_FILE, "r") as f:
    for i, line in enumerate(f):
        parts = line.strip().split()
        if len(parts) >= 2:
            fid, iid = parts[0], parts[1]
            samples.append({"FID": fid, "IID": iid})
            pid_to_indices[fid].append(i)

N = len(samples)

presence_data = []
matched_pids = []
for pid in NEW_PIDS:
    if pid in pid_to_indices:
        matched_pids.append(pid)
        presence_data.append({
            "pid": pid,
            "status": "found",
            "sample_count": len(pid_to_indices[pid]),
            "matching_field": "FID/Population",
            "notes": "Matched exactly"
        })
    else:
        presence_data.append({
            "pid": pid,
            "status": "not_found",
            "sample_count": 0,
            "matching_field": "none",
            "notes": "Not present in FAM"
        })

pd.DataFrame(presence_data).to_csv(os.path.join(OUT_DIR, "new_pid_presence_summary_260522.csv"), index=False)
print(f"Found {len(matched_pids)} out of {len(NEW_PIDS)} new PIDs.")

# 3. IID-based lookup for Region/Period
metadata = {}
try:
    with open(os.path.join(BASE_DIR, "raw/SampleInfo_260522.csv"), "r", encoding="utf-8-sig", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            iid = row.get("IID")
            if iid:
                metadata[iid] = {"region": row.get("Country", ""), "period": row.get("Period", ""), "source": "SampleInfo"}
except Exception as e:
    pass

try:
    with open(os.path.join(BASE_DIR, "raw/v66.1240K.aadr.PUB.csv"), "r", encoding="utf-8-sig", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            iid = row.get("Individual ID")
            if iid and iid not in metadata:
                group_id = row.get("Group ID", "")
                period = group_id.split("_")[-1] if "_" in group_id else group_id
                metadata[iid] = {"region": row.get("Political Entity", ""), "period": period, "source": "AADR"}
except Exception as e:
    pass

# We will just note the region/period for our samples if needed, but frequencies are calculated per-PID.

# 4. Read BIM and prepare SNP lists
known_df = pd.read_csv(KNOWN_CATALOG)
cand_df = pd.read_csv(CANDIDATES)

snp_dict = {}
for _, row in known_df.iterrows():
    snp_dict[row['rsID']] = {
        'source': 'known', 'gene': row['gene'], 'trait': row['trait_or_adaptation'], 
        'effect_allele': row['effect_or_selected_allele']
    }
for _, row in cand_df.iterrows():
    if row['rsID'] not in snp_dict:
        snp_dict[row['rsID']] = {
            'source': 'candidate', 'gene': 'gene_unknown', 'trait': 'unknown', 
            'effect_allele': 'unknown'
        }

bim_snps = []
bim_indices_to_read = []
rsid_to_idx = {}

with open(BIM_FILE, "r") as f:
    for i, line in enumerate(f):
        parts = line.strip().split()
        if len(parts) >= 6:
            rsid = parts[1]
            if rsid in snp_dict:
                bim_snps.append({"index": i, "chr": parts[0], "rsid": rsid, "pos": parts[3], "a1": parts[4], "a2": parts[5]})
                bim_indices_to_read.append(i)
                rsid_to_idx[rsid] = i

print(f"Matched {len(bim_indices_to_read)} SNPs in BIM.")

# Create the new Unified Catalog format
unified_catalog = []
for rsid, info in snp_dict.items():
    if rsid not in rsid_to_idx:
        continue
    b = bim_snps[bim_indices_to_read.index(rsid_to_idx[rsid])]
    
    label = f"{rsid}_{info['gene']}" if info['gene'] != 'gene_unknown' else rsid
    unified_catalog.append({
        "snp_label": label,
        "rsID": rsid,
        "chromosome": b["chr"],
        "position": b["pos"],
        "mapped_gene": info['gene'],
        "nearest_gene": info['gene'],
        "gene_mapping_basis": "previous_catalog" if info['source'] == 'known' else "gene_mapping_uncertain",
        "distance_to_gene_bp": "unknown",
        "variant_consequence": "unknown",
        "trait_or_adaptation": info['trait'],
        "effect_or_selected_allele": info['effect_allele'],
        "confidence": "high" if info['source'] == 'known' else "candidate",
        "notes": "Inherited from previous adaptive sweep analysis"
    })

pd.DataFrame(unified_catalog).to_csv(os.path.join(OUT_DIR, "known_selection_snp_gene_catalog_260522.csv"), index=False)

# 5. Read BED for the matched SNPs
bytes_per_snp = int(np.ceil(N / 4.0))

freq_data = []

with open(BED_FILE, "rb") as f:
    magic = f.read(3)
    
    for b_idx in bim_indices_to_read:
        b = [x for x in bim_snps if x["index"] == b_idx][0]
        rsid = b["rsid"]
        
        f.seek(3 + b_idx * bytes_per_snp)
        snp_bytes = f.read(bytes_per_snp)
        
        # Calculate frequencies per matched PID
        for pid in matched_pids:
            indices = pid_to_indices[pid]
            a2_count = 0
            called = 0
            
            for i in indices:
                byte_idx = i // 4
                bit_idx = (i % 4) * 2
                val = (snp_bytes[byte_idx] >> bit_idx) & 3
                
                if val == 0:
                    called += 1
                elif val == 2:
                    called += 1
                    a2_count += 1
                elif val == 3:
                    called += 1
                    a2_count += 2
                    
            pop_n = len(indices)
            missing = pop_n - called
            missing_rate = missing / pop_n if pop_n > 0 else 1.0
            tot_alleles = called * 2
            a2_freq = a2_count / tot_alleles if tot_alleles > 0 else "NA"
            
            freq_data.append({
                "Population": pid,
                "rsID": rsid,
                "BIM_A1": b["a1"],
                "BIM_A2": b["a2"],
                "A2_count": a2_count,
                "Total_alleles": tot_alleles,
                "A2_frequency": a2_freq,
                "Missing_rate": round(missing_rate, 4),
                "Low_Sample_Size_Warning": "Yes" if called < 5 else "No"
            })

pd.DataFrame(freq_data).to_csv(os.path.join(OUT_DIR, "new_pid_known_snp_frequencies_260522.csv"), index=False)
print("Frequencies calculated and saved.")
