import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

BASE_DIR = "/home/duri_bae/popgen_lab"
OUT_DIR = os.path.join(BASE_DIR, "output/adaptive_sweep/PID_comparison")
FIG_DIR = os.path.join(OUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

contrasts = pd.read_csv(os.path.join(OUT_DIR, "previous_targets_vs_each_new_pid_snp_contrast_260522.csv"))
catalog = pd.read_csv(os.path.join(OUT_DIR, "known_selection_snp_gene_catalog_260522.csv"))
candidates = pd.read_csv(os.path.join(OUT_DIR, "new_pid_characteristic_snp_candidates_260522.csv"))
presence = pd.read_csv(os.path.join(OUT_DIR, "new_pid_presence_summary_260522.csv"))

# Filter known
known_snps = catalog[catalog['confidence'] == 'high']['rsID'].tolist()
known_contrasts = contrasts[contrasts['rsID'].isin(known_snps)]

def clean_af(val):
    if pd.isna(val) or val == "NA": return np.nan
    return float(val)

contrasts['PID_AF_num'] = contrasts['PID_AF'].apply(clean_af)
contrasts['Target_AF_num'] = contrasts['Target_AF'].apply(clean_af)
contrasts['AF_Difference_num'] = contrasts['AF_Difference'].apply(clean_af)
known_contrasts['PID_AF_num'] = known_contrasts['PID_AF'].apply(clean_af)
known_contrasts['Target_AF_num'] = known_contrasts['Target_AF'].apply(clean_af)
known_contrasts['AF_Difference_num'] = known_contrasts['AF_Difference'].apply(clean_af)

# 1. AF heatmap of known selection SNPs across new PIDs
pivot_known_af = known_contrasts.pivot(index='rsID', columns='PID', values='PID_AF_num')
target_afs = known_contrasts.drop_duplicates(subset=['rsID']).set_index('rsID')[['Target_AF_num']]
target_afs.columns = ['Previous Targets']
combined_af = pd.concat([target_afs, pivot_known_af], axis=1)

plt.figure(figsize=(14, 8))
sns.heatmap(combined_af, annot=True, cmap="YlOrRd", vmin=0, vmax=1, fmt=".2f")
plt.title("Allele Frequency of Known Selection SNPs (Including Previous Targets)")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "1_known_snp_af_heatmap_260522.png"))
plt.close()

# 2. AF heatmap comparing previous targets and each new PID (for known SNPs)
plt.figure(figsize=(14, 8))
sns.heatmap(combined_af, annot=True, cmap="YlOrRd", vmin=0, vmax=1, fmt=".2f")
plt.title("Comparison: Previous Targets vs New PIDs (Known SNPs)")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "2_target_vs_pid_af_heatmap_260522.png"))
plt.close()

# 3. Target-vs-new-PID AF contrast heatmap
pivot_diff = known_contrasts.pivot(index='rsID', columns='PID', values='AF_Difference_num')
plt.figure(figsize=(14, 8))
sns.heatmap(pivot_diff, annot=True, cmap="coolwarm", center=0, vmin=-1, vmax=1, fmt=".2f")
plt.title("Allele Frequency Contrast (PID - Target)")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "3_target_vs_pid_contrast_heatmap_260522.png"))
plt.close()

# 4. SNP-gene-trait summary image
fig, ax = plt.subplots(figsize=(24, 12))
ax.axis('tight')
ax.axis('off')
table_data = [["SNP Label", "Gene", "Trait", "Effect Allele", "Target AF", "Max Contrast", "Confidence", "Notes"]]
for _, row in catalog[catalog['confidence'] == 'high'].iterrows():
    rsid = row['rsID']
    label = row['snp_label']
    gene = row['mapped_gene']
    trait = row['trait_or_adaptation'][:30] + "..." if len(row['trait_or_adaptation'])>30 else row['trait_or_adaptation']
    effect = row['effect_or_selected_allele']
    
    tgt_val = target_afs.loc[rsid, 'Previous Targets'] if rsid in target_afs.index else np.nan
    diffs = pivot_diff.loc[rsid] if rsid in pivot_diff.index else []
    max_diff = diffs.abs().max() if len(diffs) > 0 and not diffs.isna().all() else np.nan
    
    table_data.append([
        label, gene, trait, effect, 
        f"{tgt_val:.2f}" if pd.notna(tgt_val) else "NA", 
        f"{max_diff:.2f}" if pd.notna(max_diff) else "NA", 
        row['confidence'], "Inherited from previous analysis"
    ])
table = ax.table(cellText=table_data, loc='center', cellLoc='left')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2)
plt.title("SNP-Gene-Trait Summary Table", fontsize=16)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "snp_gene_trait_summary_new_pid_260522.png"))
plt.close()

# 5. Characteristic SNP summary image for new PIDs
fig, ax = plt.subplots(figsize=(24, 14))
ax.axis('tight')
ax.axis('off')
summary_data = [["PID", "rsID", "Gene", "Candidate Label", "PID AF", "Target AF", "Contrast Dir", "Confidence"]]
# show top 30
top_cands = candidates.head(30)
for _, row in top_cands.iterrows():
    summary_data.append([
        row['PID'], row['rsID'], row['gene'], row['Candidate_Label'], 
        str(row['PID_AF']), str(row['Target_AF']), row['Contrast_Direction'], row['Confidence']
    ])
table = ax.table(cellText=summary_data, loc='center', cellLoc='left')
table.auto_set_font_size(False)
table.set_fontsize(9)
table.scale(1, 1.5)
plt.title("Characteristic SNP Summary (Top 30 Differentiated)", fontsize=16)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "characteristic_snp_summary_new_pid_260522.png"))
plt.close()

# 6. Target-vs-PID allele-frequency scatter plot
plt.figure(figsize=(10, 8))
plt.scatter(contrasts['Target_AF_num'], contrasts['PID_AF_num'], alpha=0.3, color='blue')
plt.plot([0, 1], [0, 1], 'r--')
plt.xlabel("Previous Targets Allele Frequency")
plt.ylabel("New PIDs Allele Frequency")
plt.title("Overall AF Scatter: Targets vs PIDs")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "6_target_vs_pid_af_scatter_260522.png"))
plt.close()

# 7. Missingness rate per SNP per new PID
pivot_miss = contrasts.pivot(index='rsID', columns='PID', values='PID_Missingness')
target_miss = contrasts.drop_duplicates(subset=['rsID']).set_index('rsID')[['Target_Missingness']]
target_miss.columns = ['Previous Targets']
combined_miss = pd.concat([target_miss, pivot_miss], axis=1)

plt.figure(figsize=(14, 8))
# take a sample of 50 SNPs for readability if too large
sns.heatmap(combined_miss.head(50), cmap="Reds", vmin=0, vmax=1)
plt.title("Missingness Rate per SNP (Including Previous Targets, Sample of 50)")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "7_missingness_rate_heatmap_260522.png"))
plt.close()

# 8. Sample-size summary per new PID
plt.figure(figsize=(12, 6))
sns.barplot(data=presence, x='pid', y='sample_count', palette="viridis")
plt.xticks(rotation=45, ha='right')
plt.title("Sample Size per New PID")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "8_sample_size_summary_260522.png"))
plt.close()

# 9. High-altitude & Pigmentation SNPs by PID
interest = catalog[catalog['mapped_gene'].isin(['EPAS1', 'EGLN1', 'SLC24A5', 'OCA2', 'TYR', 'EDAR'])]['rsID'].tolist()
int_contrasts = contrasts[contrasts['rsID'].isin(interest)]
pivot_int = int_contrasts.pivot(index='rsID', columns='PID', values='PID_AF_num')

target_afs_int = int_contrasts.drop_duplicates(subset=['rsID']).set_index('rsID')[['Target_AF_num']]
target_afs_int.columns = ['Previous Targets']
combined_int = pd.concat([target_afs_int, pivot_int], axis=1)

plt.figure(figsize=(14, 6))
sns.heatmap(combined_int, annot=True, cmap="YlGnBu", vmin=0, vmax=1, fmt=".2f")
plt.title("Selected Trait SNPs AF (Including Previous Targets)")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "9_selected_traits_af_heatmap_260522.png"))
plt.close()

# 10. Shared vs distinct candidate patterns
counts = candidates.groupby(['Candidate_Label']).size().reset_index(name='Count')
plt.figure(figsize=(10, 6))
sns.barplot(data=counts, x='Candidate_Label', y='Count', palette="Set2")
plt.xticks(rotation=45, ha='right')
plt.title("Frequency of Characteristic Candidate Types")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "10_shared_vs_distinct_candidates_260522.png"))
plt.close()

print(f"Successfully generated 10 visualization files in {FIG_DIR}")
