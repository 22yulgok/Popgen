import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.filterwarnings("ignore")

BASE_DIR = "/home/duri_bae/popgen_lab"
OUT_DIR = os.path.join(BASE_DIR, "output/adaptive_sweep/selection_tests")
FIG_DIR = os.path.join(OUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

freq_df = pd.read_csv(os.path.join(OUT_DIR, "selected_snp_population_frequencies_260522.csv"))
grp_summary = pd.read_csv(os.path.join(OUT_DIR, "selection_test_group_summary_260522.csv"))
pos_df = pd.read_csv(os.path.join(OUT_DIR, "positive_selection_frequency_tests_260522.csv"))
fst_pbs = pd.read_csv(os.path.join(OUT_DIR, "selected_snp_fst_pbs_summary_260522.csv"))
window_df = pd.read_csv(os.path.join(OUT_DIR, "local_window_selection_summary_260522.csv"))
interp_df = pd.read_csv(os.path.join(OUT_DIR, "selection_test_interpretation_summary_260522.csv"))

def safe_float(x):
    try: return float(x)
    except: return np.nan
freq_df['A2_frequency_num'] = freq_df['A2_frequency'].apply(safe_float)

pop_class = dict(zip(grp_summary['population'], grp_summary['group_class']))
freq_df['Group_Class'] = freq_df['Group'].map(pop_class)

# 1. Selected SNP AF heatmap across all comparison groups (sampled top SNPs)
pivot_af = freq_df.pivot(index='rsID', columns='Group', values='A2_frequency_num')
pivot_af = pivot_af.dropna(thresh=20) # Keep SNPs typed in most groups
plt.figure(figsize=(12, 10))
sns.heatmap(pivot_af.head(40), cmap="YlGnBu", vmin=0, vmax=1)
plt.title("Allele Frequency Across All Groups")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "1_selected_snp_af_heatmap_all_groups_260522.png"))
plt.close()

# 2. Allele-frequency heatmap by group class
class_af = freq_df.groupby(['rsID', 'Group_Class'])['A2_frequency_num'].mean().reset_index()
pivot_class = class_af.pivot(index='rsID', columns='Group_Class', values='A2_frequency_num').dropna()
plt.figure(figsize=(8, 10))
sns.heatmap(pivot_class.head(40), annot=True, cmap="YlGnBu", vmin=0, vmax=1)
plt.title("Average Allele Frequency by Group Class")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "2_allele_frequency_heatmap_by_group_class_260522.png"))
plt.close()

# 3. Basal vs European vs Eastern Eurasian Boxplot
plt.figure(figsize=(10, 6))
sns.boxplot(data=freq_df.dropna(subset=['A2_frequency_num']), x='Group_Class', y='A2_frequency_num')
plt.title("Global AF Distribution by Broad Ancestry Class")
plt.xticks(rotation=15)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "3_basal_vs_european_vs_eastern_comparison_260522.png"))
plt.close()

# 4. Pairwise allele-frequency contrast heatmap (Max Delta AF between pop pairs)
if not pos_df.empty:
    pair_max = pos_df.groupby(['population_1', 'population_2'])['frequency_difference'].apply(lambda x: x.abs().max()).reset_index()
    pair_pivot = pair_max.pivot(index='population_1', columns='population_2', values='frequency_difference')
    plt.figure(figsize=(12, 12))
    sns.heatmap(pair_pivot, cmap="Reds", vmin=0, vmax=1)
    plt.title("Maximum AF Differentiation Between Population Pairs")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "4_pairwise_allele_frequency_contrast_heatmap_260522.png"))
    plt.close()

# 5. Positive-selection-compatible SNP summary
if not pos_df.empty:
    top_pos = pos_df[pos_df['positive_selection_compatibility_label'] == 'positive_selection_compatible']
    if not top_pos.empty:
        idx = top_pos.groupby('rsID')['frequency_difference'].apply(lambda x: x.abs().idxmax())
        top_snps = top_pos.loc[idx].sort_values(by='frequency_difference', ascending=False, key=abs).head(15)
        top_snps['Pair'] = top_snps['population_1'] + " vs " + top_snps['population_2']
        
        plt.figure(figsize=(10,6))
        sns.barplot(data=top_snps, x='rsID', y='frequency_difference', hue='Pair', dodge=False)
        plt.xticks(rotation=45, ha='right')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.title("Top Positive Selection Compatible SNPs (Max Pairwise Contrast)")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "5_positive_selection_compatible_summary_260522.png"))
        plt.close()

# 6. FST or PBS summary plot
if not fst_pbs.empty:
    plt.figure(figsize=(8,6))
    sns.scatterplot(data=fst_pbs, x='FST_Pop_European', y='PBS_Score', hue='Population', legend=False, alpha=0.5)
    plt.xlabel("FST (Population vs European)")
    plt.ylabel("PBS Score")
    plt.title("PBS vs FST Distribution Across All Target & PID Groups")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "6_fst_pbs_summary_plot_260522.png"))
    plt.close()

# 7. Local window outlier plot
if not window_df.empty:
    top_w = window_df.sort_values(by="core_Pop_vs_European_FST", ascending=False).head(20)
    plt.figure(figsize=(12, 6))
    bar_width = 0.35
    x = np.arange(len(top_w))
    plt.bar(x, top_w['core_Pop_vs_European_FST'], bar_width, label='Core SNP FST')
    plt.bar(x + bar_width, top_w['window_mean_FST'], bar_width, label='Window Mean FST')
    labels = [f"{r} ({p})" for r, p in zip(top_w['rsID'], top_w['population_tested'])]
    plt.xticks(x + bar_width/2, labels, rotation=45, ha='right')
    plt.legend()
    plt.title("Local FST Outliers (Core vs ±500kb Mean)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "7_local_window_outlier_summary_260522.png"))
    plt.close()

# 8. Temporal AF trend mock
plt.figure(figsize=(8, 6))
plt.text(0.5, 0.5, "Insufficient longitudinal IID dates to plot robust Temporal AF.\n(See temporal test caveats)", ha='center', va='center', fontsize=12)
plt.axis('off')
plt.title("Temporal AF Trend")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "8_temporal_af_trend_mock_260522.png"))
plt.close()

# 9. P-value/FDR summary plot
if not pos_df.empty:
    plt.figure(figsize=(8,6))
    plt.hist(pos_df['p_value'].dropna(), bins=50, color='skyblue', edgecolor='black')
    plt.axvline(0.05, color='r', linestyle='--')
    plt.yscale('log')
    plt.xlabel("P-value")
    plt.ylabel("Log Frequency")
    plt.title("Fisher Exact Test P-value Distribution (All Pairs)")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "9_pvalue_fdr_summary_plot_260522.png"))
    plt.close()

# 10. Missingness and sample-size caution
plt.figure(figsize=(8,6))
plt.scatter(freq_df['Sample_Size'], freq_df['Missing_rate'], alpha=0.3, color='orange')
plt.xlabel("Total Sample Size")
plt.ylabel("Missing Rate")
plt.title("Missingness vs Sample Size Across Groups")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "10_missingness_sample_size_caution_260522.png"))
plt.close()

# 11. SNP-gene-selection interpretation summary image
fig, ax = plt.subplots(figsize=(24, 10))
ax.axis('tight')
ax.axis('off')

# Compile display table
table_data = [["rsID", "Gene", "Trait", "Allele", "Highest Freq Pop", "Strongest Contrast Pair", "Overall Label", "Outlier?", "Caveats"]]
for _, row in interp_df.head(25).iterrows():
    trait = str(row['trait_or_adaptation'])[:25] + "..." if len(str(row['trait_or_adaptation']))>25 else str(row['trait_or_adaptation'])
    contrast = str(row['strongest_population_contrast'])[:35] + "..." if len(str(row['strongest_population_contrast']))>35 else str(row['strongest_population_contrast'])
    caveats = str(row['caveats'])[:20] + "..." if len(str(row['caveats']))>20 else str(row['caveats'])
    table_data.append([
        row['rsID'], str(row['mapped_gene']), trait, str(row['tested_allele']),
        str(row['highest_frequency_population']), contrast, str(row['overall_selection_compatibility']), 
        str(row['is_local_window_outlier']), caveats
    ])

table = ax.table(cellText=table_data, loc='center', cellLoc='left')
table.auto_set_font_size(False)
table.set_fontsize(8)
table.scale(1, 1.5)
plt.title("SNP-Gene Positive Selection Test Summary", fontsize=16)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "snp_gene_positive_selection_test_summary_260522.png"))
plt.close()

print("Generated 11 visualization plots in output directory.")
