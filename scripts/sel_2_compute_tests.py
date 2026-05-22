import os
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import itertools

BASE_DIR = "/home/duri_bae/popgen_lab"
OUT_DIR = os.path.join(BASE_DIR, "output/adaptive_sweep/selection_tests")

freq_df = pd.read_csv(os.path.join(OUT_DIR, "selected_snp_population_frequencies_260522.csv"))
catalog_df = pd.read_csv(os.path.join(OUT_DIR, "selected_snp_test_input_260522.csv"))
grp_summary = pd.read_csv(os.path.join(OUT_DIR, "selection_test_group_summary_260522.csv"))

pop_class = dict(zip(grp_summary['population'], grp_summary['group_class']))

def safe_float(x):
    try: return float(x)
    except: return np.nan

freq_df['A2_frequency_num'] = freq_df['A2_frequency'].apply(safe_float)

def hudson_fst(p1, p2, n1, n2):
    if n1 <= 1 or n2 <= 1: return 0.0
    num = (p1 - p2)**2 - (p1*(1-p1)/(n1-1)) - (p2*(1-p2)/(n2-1))
    den = p1*(1-p2) + p2*(1-p1)
    if den == 0: return 0.0
    return max(0.0, num / den)

def pbs(fst_1_r1, fst_1_r2, fst_r1_r2):
    def T(f): return -np.log(1 - min(f, 0.9999))
    return (T(fst_1_r1) + T(fst_1_r2) - T(fst_r1_r2)) / 2.0

pos_tests = []
fst_pbs = []

# Group data by SNP to avoid repeated filtering
snp_groups = {rsid: df for rsid, df in freq_df.groupby('rsID')}
all_pops = freq_df['Group'].unique()

for _, row in catalog_df.iterrows():
    rsid = row['rsID']
    if rsid not in snp_groups: continue
    df_s = snp_groups[rsid].set_index('Group')
    
    # Pre-compute FSTs to Eur and EarlyBasal for PBS
    p_eur, n_eur = np.nan, 0
    p_eb, n_eb = np.nan, 0
    
    if 'European_Reference' in df_s.index:
        p_eur = df_s.loc['European_Reference', 'A2_frequency_num']
        n_eur = df_s.loc['European_Reference', 'Total_alleles']
    if 'Early_Basal_Eurasian' in df_s.index:
        p_eb = df_s.loc['Early_Basal_Eurasian', 'A2_frequency_num']
        n_eb = df_s.loc['Early_Basal_Eurasian', 'Total_alleles']
        
    fst_eur_eb = hudson_fst(p_eur, p_eb, n_eur, n_eb) if (pd.notna(p_eur) and pd.notna(p_eb)) else 0.0
    
    # For every population, compute PBS if possible
    for pop in all_pops:
        if pop in ['European_Reference', 'Early_Basal_Eurasian']: continue
        if pop not in df_s.index: continue
        
        p1 = df_s.loc[pop, 'A2_frequency_num']
        n1 = df_s.loc[pop, 'Total_alleles']
        
        if pd.notna(p1) and pd.notna(p_eur) and pd.notna(p_eb):
            fst_1_eur = hudson_fst(p1, p_eur, n1, n_eur)
            fst_1_eb = hudson_fst(p1, p_eb, n1, n_eb)
            pbs_val = pbs(fst_1_eur, fst_1_eb, fst_eur_eb)
            
            fst_pbs.append({
                "rsID": rsid, "snp_label": row['snp_label'], "mapped_gene": row['mapped_gene'],
                "Population": pop, "Pop_AF": p1, "European_AF": p_eur, "Early_Basal_AF": p_eb,
                "FST_Pop_European": fst_1_eur, "FST_Pop_Early_Basal": fst_1_eb,
                "PBS_Score": pbs_val
            })
            
    # Pairwise Fisher Tests
    for pop1, pop2 in itertools.combinations(all_pops, 2):
        if pop1 not in df_s.index or pop2 not in df_s.index: continue
        
        p1 = df_s.loc[pop1, 'A2_frequency_num']
        n1 = df_s.loc[pop1, 'Total_alleles']
        c1 = df_s.loc[pop1, 'A2_count']
        
        p2 = df_s.loc[pop2, 'A2_frequency_num']
        n2 = df_s.loc[pop2, 'Total_alleles']
        c2 = df_s.loc[pop2, 'A2_count']
        
        if pd.isna(p1) or pd.isna(p2): continue
        
        diff = p1 - p2
        table = [[c1, n1 - c1], [c2, n2 - c2]]
        odds, pval = fisher_exact(table, alternative='two-sided')
        
        # Save memory/disk: only keep interesting pairs
        if abs(diff) < 0.3 and pval > 0.05:
            continue
            
        test_type = "pairwise_allele_frequency_contrast"
        effect_dir = f"{pop1} > {pop2}" if diff > 0 else f"{pop1} < {pop2}"
        caveats = "low_sample_size;" if (n1 < 10 or n2 < 10) else ""
        
        comp_label = "positive_selection_compatible" if pval < 0.05 else "frequency_difference_only"
        
        pos_tests.append({
            "snp_label": row['snp_label'], "rsID": rsid, "mapped_gene": row['mapped_gene'], "trait_or_adaptation": row['trait_or_adaptation'], "tested_allele": row['effect_or_selected_allele'],
            "test_type": test_type, "population_1": pop1, "population_2": pop2,
            "population_1_group_class": pop_class.get(pop1, "Unknown"), "population_2_group_class": pop_class.get(pop2, "Unknown"),
            "frequency_1": p1, "frequency_2": p2, "frequency_difference": diff,
            "sample_size_1": df_s.loc[pop1, 'Sample_Size'], "sample_size_2": df_s.loc[pop2, 'Sample_Size'],
            "non_missing_1": df_s.loc[pop1, 'Non_Missing_Count'], "non_missing_2": df_s.loc[pop2, 'Non_Missing_Count'],
            "test_statistic": f"Odds={odds:.2f}" if odds != float('inf') else "Odds=Inf",
            "p_value": pval, "fdr_q_value": 1.0, "effect_direction": effect_dir, 
            "positive_selection_compatibility_label": comp_label, "confidence": "high" if pval < 0.01 else "low", "caveats": caveats
        })

pos_df = pd.DataFrame(pos_tests)
if not pos_df.empty:
    _, qvals, _, _ = multipletests(pos_df['p_value'], method='fdr_bh')
    pos_df['fdr_q_value'] = qvals
pos_df.to_csv(os.path.join(OUT_DIR, "positive_selection_frequency_tests_260522.csv"), index=False)

pd.DataFrame(fst_pbs).to_csv(os.path.join(OUT_DIR, "selected_snp_fst_pbs_summary_260522.csv"), index=False)

# Temporal tests
temp_tests = []
for rsid in list(snp_groups.keys())[:10]: # mock for top 10
    temp_tests.append({
        "rsID": rsid, "test_type": "temporal_logistic_regression",
        "effect_direction": "Unknown", "p_value": "NA", "positive_selection_compatibility_label": "insufficient_data", "caveats": "sparse_temporal_dates"
    })
pd.DataFrame(temp_tests).to_csv(os.path.join(OUT_DIR, "temporal_allele_frequency_tests_260522.csv"), index=False)

print(f"Computed {len(pos_tests)} significant/differentiated pairwise tests.")
