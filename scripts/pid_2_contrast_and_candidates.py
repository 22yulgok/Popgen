import os
import pandas as pd
import numpy as np

BASE_DIR = "/home/duri_bae/popgen_lab"
OUT_DIR = os.path.join(BASE_DIR, "output/adaptive_sweep/PID_comparison")
REPORT_FILE = os.path.join(BASE_DIR, "output/adaptive_sweep/adaptive_sweep_report_260522.md")

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

# Load data
old_freqs = pd.read_csv(os.path.join(BASE_DIR, "output/adaptive_sweep/known_selection_snp_population_frequencies_260522.csv"))
new_freqs = pd.read_csv(os.path.join(OUT_DIR, "new_pid_known_snp_frequencies_260522.csv"))
catalog = pd.read_csv(os.path.join(OUT_DIR, "known_selection_snp_gene_catalog_260522.csv"))

# Calculate pooled Target AF for each SNP
target_df = old_freqs[old_freqs['Population'].isin(targets)]
pooled_target = target_df.groupby('rsID').agg({'A2_count': 'sum', 'Total_alleles': 'sum', 'Missing_rate': 'mean'})
pooled_target['pooled_target_AF'] = np.where(pooled_target['Total_alleles'] > 0, pooled_target['A2_count'] / pooled_target['Total_alleles'], np.nan)
pooled_target = pooled_target.reset_index()

# Contrast each new PID against the pooled target AF
contrasts = []
candidates = []

for _, row in new_freqs.iterrows():
    rsid = row['rsID']
    pid = row['Population']
    pid_af = pd.to_numeric(row['A2_frequency'], errors='coerce')
    pid_miss = row['Missing_rate']
    pid_n = row['Total_alleles'] / 2 if pd.notna(row['Total_alleles']) else 0
    
    tgt_row = pooled_target[pooled_target['rsID'] == rsid]
    if tgt_row.empty:
        tgt_af = np.nan
        tgt_miss = np.nan
        tgt_n = 0
    else:
        tgt_af = tgt_row['pooled_target_AF'].values[0]
        tgt_miss = tgt_row['Missing_rate'].values[0]
        tgt_n = tgt_row['Total_alleles'].values[0] / 2
        
    diff = pid_af - tgt_af if (pd.notna(pid_af) and pd.notna(tgt_af)) else np.nan
    direction = "PID > Target" if diff > 0 else "Target > PID" if diff < 0 else "Equal/Unknown"
    
    # Label logic
    label = "not_candidate"
    if pid_miss > 0.5 or pd.isna(pid_af) or pid_n < 2:
        label = "insufficient_data"
        confidence = "low"
    else:
        confidence = "high" if pid_n >= 5 else "caution_low_sample_size"
        taf = tgt_af if pd.notna(tgt_af) else 0.0
        
        if pid_af >= 0.8 and taf <= 0.2:
            label = "pid_high_target_low_candidate"
        elif taf >= 0.8 and pid_af <= 0.2:
            label = "target_high_pid_low_candidate"
        elif pid_af >= 0.8 and taf >= 0.8:
            label = "shared_high_frequency_candidate"
        elif pid_af >= 0.9:
            label = "pid_fixed_candidate"
        elif pid_af >= 0.8:
            label = "pid_high_frequency_candidate"
            
    # Lookup catalog info
    cat_row = catalog[catalog['rsID'] == rsid]
    gene = cat_row['mapped_gene'].values[0] if not cat_row.empty else "gene_unknown"
    trait = cat_row['trait_or_adaptation'].values[0] if not cat_row.empty else "unknown"
    
    contrasts.append({
        "PID": pid,
        "rsID": rsid,
        "gene": gene,
        "trait": trait,
        "PID_AF": round(pid_af, 4) if pd.notna(pid_af) else "NA",
        "Target_AF": round(tgt_af, 4) if pd.notna(tgt_af) else "NA",
        "AF_Difference": round(diff, 4) if pd.notna(diff) else "NA",
        "Contrast_Direction": direction,
        "PID_N": pid_n,
        "Target_N": tgt_n,
        "PID_Missingness": round(pid_miss, 4),
        "Target_Missingness": round(tgt_miss, 4),
        "Candidate_Label": label,
        "Confidence": confidence
    })
    
    if label not in ["not_candidate", "insufficient_data"]:
        candidates.append(contrasts[-1])

pd.DataFrame(contrasts).to_csv(os.path.join(OUT_DIR, "previous_targets_vs_each_new_pid_snp_contrast_260522.csv"), index=False)
pd.DataFrame(candidates).to_csv(os.path.join(OUT_DIR, "new_pid_characteristic_snp_candidates_260522.csv"), index=False)

print(f"Computed contrasts for {len(new_freqs)} SNP-PID pairs.")
print(f"Identified {len(candidates)} characteristic candidate occurrences across the new PIDs.")
