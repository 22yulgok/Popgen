import os
import csv
import re
import numpy as np
import collections
import pandas as pd

BASE_DIR = "/home/duri_bae/popgen_lab"
OUT_DIR = os.path.join(BASE_DIR, "output/adaptive_sweep/selection_tests")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, "figures"), exist_ok=True)

REPORT_FILE = os.path.join(BASE_DIR, "output/adaptive_sweep/adaptive_sweep_report_260522.md")
KNOWN_CATALOG = os.path.join(BASE_DIR, "output/adaptive_sweep/PID_comparison/known_selection_snp_gene_catalog_260522.csv")
CAND_CATALOG = os.path.join(BASE_DIR, "output/adaptive_sweep/PID_comparison/new_pid_characteristic_snp_candidates_260522.csv")

FAM_FILE = os.path.join(BASE_DIR, "raw/Tibetan.fstat.1240K.260508.fam")
BIM_FILE = os.path.join(BASE_DIR, "raw/Tibetan.fstat.1240K.260508.bim")
BED_FILE = os.path.join(BASE_DIR, "raw/Tibetan.fstat.1240K.260508.bed")

NEW_PIDS = ['Haneyi', 'Xingyi_EN', 'BGD_EN', 'YR_LN', 'Upper_YR_LN', 'Shandong_EN', 'Xiaojingshan', 'Yangshao_UYR', 'Shimao_group', 'STM_EN', 'Boisman_MN', 'Yumin', 'HMMH_MN', 'DevilsCave_N', 'Qihe3', 'Liangdao', 'Xitoucun', 'Xueshan']
EARLY_BASAL = ['Ust_Ishim', 'BachoKiro_IUP', 'BachoKiro_UP', 'ZlatyKun_IUP.SG']
EUROPEAN = ['Sunghir']
EASTERN_BASAL = ['Tianyuan', 'AR33K', 'Onge.DG', 'McColl_SEA_GR1']

# 1. Infer previous targets
targets = set()
with open(REPORT_FILE, 'r') as f:
    for line in f:
        if line.startswith("- **Target Populations Found**:"):
            pops_str = line.split(":", 1)[1].strip()
            for p in pops_str.split(","):
                p = p.strip()
                if p: targets.add(p)

print(f"Inferred {len(targets)} independent previous targets.")

# Combine SNPs
catalog_df = pd.read_csv(KNOWN_CATALOG)
cand_df = pd.read_csv(CAND_CATALOG)

test_input = []
seen_rsid = set()

# Known SNPs
for _, row in catalog_df.iterrows():
    rsid = row['rsID']
    if rsid not in seen_rsid:
        seen_rsid.add(rsid)
        test_input.append({
            "snp_label": row['snp_label'],
            "rsID": rsid,
            "chromosome": row['chromosome'],
            "position": row['position'],
            "mapped_gene": row['mapped_gene'],
            "trait_or_adaptation": row['trait_or_adaptation'],
            "effect_or_selected_allele": row['effect_or_selected_allele'],
            "ancestral_allele": "unknown",
            "derived_allele": "unknown",
            "previous_category": "known_selection",
            "source_previous_file": "known_selection_snp_gene_catalog_260522.csv",
            "source_section_or_table": "catalog",
            "confidence": row['confidence'],
            "notes": row['notes']
        })

# Candidate SNPs
for _, row in cand_df.iterrows():
    rsid = row['rsID']
    if rsid not in seen_rsid:
        seen_rsid.add(rsid)
        test_input.append({
            "snp_label": row.get('snp_label', rsid),
            "rsID": rsid,
            "chromosome": "unknown", # will get from bim
            "position": "unknown", # will get from bim
            "mapped_gene": row.get('gene', 'unknown'),
            "trait_or_adaptation": "unknown",
            "effect_or_selected_allele": "unknown",
            "ancestral_allele": "unknown",
            "derived_allele": "unknown",
            "previous_category": row.get('Candidate_Label', 'differentiated_candidate'),
            "source_previous_file": "new_pid_characteristic_snp_candidates_260522.csv",
            "source_section_or_table": "candidates",
            "confidence": row.get('Confidence', 'low'),
            "notes": "Inherited from PID comparison candidates"
        })

# Assign Groups
samples = []
group_to_indices = collections.defaultdict(list)
pop_group_class = {}

with open(FAM_FILE, "r") as f:
    for i, line in enumerate(f):
        parts = line.strip().split()
        if len(parts) >= 2:
            pop, iid = parts[0], parts[1]
            samples.append({"FID": pop, "IID": iid})
            
            if pop in targets:
                group_to_indices[pop].append(i)
                pop_group_class[pop] = "Previous_Target"
            elif pop in NEW_PIDS:
                group_to_indices[pop].append(i)
                pop_group_class[pop] = "Eastern_Eurasian_PID"
            elif pop in EARLY_BASAL:
                group_to_indices["Early_Basal_Eurasian"].append(i)
            elif pop in EUROPEAN:
                group_to_indices["European_Reference"].append(i)
            elif pop in EASTERN_BASAL:
                group_to_indices["Eastern_Basal_Eurasian"].append(i)

pop_group_class["Early_Basal_Eurasian"] = "Deep_Reference"
pop_group_class["European_Reference"] = "Deep_Reference"
pop_group_class["Eastern_Basal_Eurasian"] = "Deep_Reference"

grp_summary = []
for g, idxs in group_to_indices.items():
    grp_summary.append({
        "population": g,
        "group_class": pop_group_class[g],
        "sample_count": len(idxs),
        "selection_reason": "Inferred from previous run or explicitly defined baseline",
        "status": "Found",
        "notes": ""
    })
pd.DataFrame(grp_summary).to_csv(os.path.join(OUT_DIR, "selection_test_group_summary_260522.csv"), index=False)

# Match SNPs
rsid_to_idx = {}
bim_snps = []
bim_indices_to_read = []

test_rsids = set(seen_rsid)

with open(BIM_FILE, "r") as f:
    for i, line in enumerate(f):
        parts = line.strip().split()
        if len(parts) >= 6:
            rsid = parts[1]
            if rsid in test_rsids:
                bim_snps.append({"index": i, "chr": parts[0], "rsid": rsid, "pos": parts[3], "a1": parts[4], "a2": parts[5]})
                bim_indices_to_read.append(i)
                rsid_to_idx[rsid] = {"index": i, "chr": parts[0], "pos": parts[3]}

match_status = []
for ti in test_input:
    r = ti['rsID']
    if r in rsid_to_idx:
        match_status.append({"rsID": r, "matched": "Yes", "bim_index": rsid_to_idx[r]['index']})
        if ti['chromosome'] == 'unknown':
            ti['chromosome'] = rsid_to_idx[r]['chr']
            ti['position'] = rsid_to_idx[r]['pos']
    else:
        match_status.append({"rsID": r, "matched": "No", "bim_index": "NA"})

pd.DataFrame(test_input).to_csv(os.path.join(OUT_DIR, "selected_snp_test_input_260522.csv"), index=False)
pd.DataFrame(match_status).to_csv(os.path.join(OUT_DIR, "selected_snp_match_status_260522.csv"), index=False)

N = len(samples)
bytes_per_snp = int(np.ceil(N / 4.0))
freq_data = []

with open(BED_FILE, "rb") as f:
    magic = f.read(3)
    for b in bim_snps:
        b_idx = b["index"]
        rsid = b["rsid"]
        
        f.seek(3 + b_idx * bytes_per_snp)
        snp_bytes = f.read(bytes_per_snp)
        
        for g, indices in group_to_indices.items():
            a2_count = 0
            called = 0
            
            for i in indices:
                val = (snp_bytes[i // 4] >> ((i % 4) * 2)) & 3
                if val == 0: called += 1
                elif val == 2: called += 1; a2_count += 1
                elif val == 3: called += 1; a2_count += 2
                    
            pop_n = len(indices)
            missing = pop_n - called
            missing_rate = missing / pop_n if pop_n > 0 else 1.0
            tot_alleles = called * 2
            a2_freq = a2_count / tot_alleles if tot_alleles > 0 else "NA"
            
            freq_data.append({
                "Group": g,
                "rsID": rsid,
                "A2_count": a2_count,
                "Total_alleles": tot_alleles,
                "A2_frequency": a2_freq,
                "Missing_rate": round(missing_rate, 4),
                "Sample_Size": pop_n,
                "Non_Missing_Count": called
            })

pd.DataFrame(freq_data).to_csv(os.path.join(OUT_DIR, "selected_snp_population_frequencies_260522.csv"), index=False)
print(f"Frequencies calculated for {len(group_to_indices)} broad groups.")

# Dates
metadata_dates = {}
try:
    with open(os.path.join(BASE_DIR, "raw/v66.1240K.aadr.PUB.csv"), "r", encoding="utf-8-sig", errors="ignore") as f:
        reader = csv.DictReader(f)
        date_col = next((c for c in reader.fieldnames if c and "Date mean in BP" in c), None)
        
        for row in reader:
            iid = row.get("Individual ID")
            date_str = row.get(date_col) if date_col else None
            if iid and date_str:
                try: metadata_dates[iid] = -float(date_str)
                except: pass
except Exception: pass

date_out = [{"IID": i, "Numeric_Date": d} for i, d in metadata_dates.items()]
pd.DataFrame(date_out).to_csv(os.path.join(OUT_DIR, "iid_dates_260522.csv"), index=False)
print(f"Extracted dates for {len(date_out)} individuals.")
