import os
import pandas as pd
import numpy as np
import collections

BASE_DIR = "/home/duri_bae/popgen_lab"
OUT_DIR = os.path.join(BASE_DIR, "output/adaptive_sweep/selection_tests")

# We will analyze the top 20 selected SNPs based on the strongest Pop vs European FST.
fst_df = pd.read_csv(os.path.join(OUT_DIR, "selected_snp_fst_pbs_summary_260522.csv"))
# Get the max FST vs European per SNP
idx = fst_df.groupby('rsID')['FST_Pop_European'].idxmax()
top_snps = fst_df.loc[idx].sort_values(by="FST_Pop_European", ascending=False).head(20)

top_rsids = set(top_snps['rsID'].tolist())
rsid_to_pop = dict(zip(top_snps['rsID'], top_snps['Population']))

BIM_FILE = os.path.join(BASE_DIR, "raw/Tibetan.fstat.1240K.260508.bim")
BED_FILE = os.path.join(BASE_DIR, "raw/Tibetan.fstat.1240K.260508.bed")
FAM_FILE = os.path.join(BASE_DIR, "raw/Tibetan.fstat.1240K.260508.fam")

# Get Group indices for all pops in rsid_to_pop, plus European_Reference
needed_pops = set(rsid_to_pop.values())
needed_pops.add("European_Reference")

group_to_indices = collections.defaultdict(list)

N = 0
with open(FAM_FILE, "r") as f:
    for i, line in enumerate(f):
        pop = line.strip().split()[0]
        if pop in needed_pops:
            group_to_indices[pop].append(i)
        N += 1

bytes_per_snp = int(np.ceil(N / 4.0))

snps_by_chr = collections.defaultdict(list)
core_snps = {}

with open(BIM_FILE, "r") as f:
    for i, line in enumerate(f):
        parts = line.strip().split()
        if len(parts) >= 4:
            chrom, rsid, pos = parts[0], parts[1], int(parts[3])
            snps_by_chr[chrom].append({"index": i, "rsid": rsid, "pos": pos})
            if rsid in top_rsids:
                core_snps[rsid] = {"index": i, "chr": chrom, "pos": pos}

def hudson_fst(p1, p2, n1, n2):
    if n1 <= 1 or n2 <= 1: return 0.0
    num = (p1 - p2)**2 - (p1*(1-p1)/(n1-1)) - (p2*(1-p2)/(n2-1))
    den = p1*(1-p2) + p2*(1-p1)
    if den == 0: return 0.0
    return max(0.0, num / den)

WINDOW_SIZE = 500000

window_summaries = []
outlier_status_map = {}

with open(BED_FILE, "rb") as bed:
    for rsid, core in core_snps.items():
        chrom = core["chr"]
        core_pos = core["pos"]
        pop = rsid_to_pop[rsid]
        
        window_snps = [s for s in snps_by_chr[chrom] if abs(s["pos"] - core_pos) <= WINDOW_SIZE]
        if not window_snps: continue
        
        local_fsts = []
        core_fst = 0.0
        
        for w_snp in window_snps:
            b_idx = w_snp["index"]
            bed.seek(3 + b_idx * bytes_per_snp)
            snp_bytes = bed.read(bytes_per_snp)
            
            p_t, p_e = 0, 0
            n_t, n_e = 0, 0
            
            for g in [pop, "European_Reference"]:
                a2_count = 0
                called = 0
                for i in group_to_indices[g]:
                    val = (snp_bytes[i // 4] >> ((i % 4) * 2)) & 3
                    if val == 0: called += 1
                    elif val == 2: called += 1; a2_count += 1
                    elif val == 3: called += 1; a2_count += 2
                if g == pop:
                    p_t = a2_count / (called*2) if called > 0 else 0
                    n_t = called*2
                else:
                    p_e = a2_count / (called*2) if called > 0 else 0
                    n_e = called*2
            
            fst = hudson_fst(p_t, p_e, n_t, n_e)
            local_fsts.append(fst)
            if w_snp["rsid"] == rsid:
                core_fst = fst
                
        if len(local_fsts) > 10:
            percentile = (sum(1 for x in local_fsts if x <= core_fst) / len(local_fsts)) * 100
        else:
            percentile = "NA"
            
        outlier = "Yes" if (percentile != "NA" and percentile > 95) else "No"
        outlier_status_map[rsid] = outlier
        
        window_summaries.append({
            "rsID": rsid,
            "chr": chrom,
            "core_pos": core_pos,
            "population_tested": pop,
            "snps_in_window": len(local_fsts),
            "core_Pop_vs_European_FST": core_fst,
            "window_mean_FST": np.mean(local_fsts),
            "percentile_in_window": percentile,
            "is_outlier_95th": outlier
        })

pd.DataFrame(window_summaries).to_csv(os.path.join(OUT_DIR, "local_window_selection_summary_260522.csv"), index=False)

# Compile overall interpretation summary for ALL SNPs (taking the strongest contrast logic)
catalog_df = pd.read_csv(os.path.join(OUT_DIR, "selected_snp_test_input_260522.csv"))
interp = []

pos_tests = pd.read_csv(os.path.join(OUT_DIR, "positive_selection_frequency_tests_260522.csv"))
pos_grouped = pos_tests.groupby('rsID')

for _, row in catalog_df.iterrows():
    rsid = row['rsID']
    outlier = outlier_status_map.get(rsid, "Not Tested")
    
    lbl = "frequency_difference_only"
    caveats = str(row['notes']) + "; unphased"
    pop_highest = "Unknown"
    strongest_contrast = "Unknown"
    
    if rsid in pos_grouped.groups:
        grp = pos_grouped.get_group(rsid)
        if any(grp['positive_selection_compatibility_label'] == 'positive_selection_compatible'):
            lbl = "positive_selection_compatible"
            if outlier == "Yes": lbl += "_with_local_support"
            
        # Get highest freq pop
        try:
            freq_df = pd.read_csv(os.path.join(OUT_DIR, "selected_snp_population_frequencies_260522.csv"))
            f_rsid = freq_df[freq_df['rsID'] == rsid]
            if not f_rsid.empty:
                pop_highest = f_rsid.loc[f_rsid['A2_frequency'].astype(float).idxmax()]['Group']
        except:
            pass
            
        # Get strongest contrast
        strongest = grp.loc[grp['frequency_difference'].abs().idxmax()]
        strongest_contrast = f"{strongest['population_1']} vs {strongest['population_2']} (Delta={strongest['frequency_difference']:.2f})"
        
    else:
        lbl = "insufficient_data"
        caveats += "; not matched or sparse"

    interp.append({
        "rsID": rsid,
        "snp_label": row['snp_label'],
        "mapped_gene": row['mapped_gene'],
        "trait_or_adaptation": row['trait_or_adaptation'],
        "tested_allele": row['effect_or_selected_allele'],
        "highest_frequency_population": pop_highest,
        "strongest_population_contrast": strongest_contrast,
        "overall_selection_compatibility": lbl,
        "is_local_window_outlier": outlier,
        "caveats": caveats
    })

pd.DataFrame(interp).to_csv(os.path.join(OUT_DIR, "selection_test_interpretation_summary_260522.csv"), index=False)
print("Local window testing and interpretation generation completed.")
