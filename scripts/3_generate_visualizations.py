import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import collections

OUTPUT_DIR = "/home/duri_bae/popgen_lab/output/adaptive_sweep"
FIG_DIR = os.path.join(OUTPUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

TARGET_POPS = {'GBSL', 'GBSL_old1', 'GBSL_old2', 'Yushu2.8k', 'Zongri5.1k', 'Zongri4.7k', 'Zongri4.1k', 'Zongri4.5k', 'Shannan3k', 'Chamdo2.8k_1', 'Mbc3.5k_G0', 'Mbc3.5k_outlier', 'Mbc4.4k_G0', 'Mbc4.4k_G1', 'Mbc4.4k_G2', 'Mbc4k_G2', 'SUI001.WGS'}

# Load data
match_status = pd.read_csv(os.path.join(OUTPUT_DIR, "known_selection_snp_match_status_260522.csv"))
freqs = pd.read_csv(os.path.join(OUTPUT_DIR, "known_selection_snp_population_frequencies_260522.csv"))
candidates = pd.read_csv(os.path.join(OUTPUT_DIR, "high_frequency_polymorphic_candidates_260522.csv"))
catalog = pd.read_csv(os.path.join(OUTPUT_DIR, "known_selection_snp_catalog_260522.csv"))

# Merge catalog info into freqs
freqs = freqs.merge(catalog[['rsID', 'gene', 'trait_or_adaptation']], on='rsID', how='left')

# Convert AF to numeric
freqs['A2_frequency'] = pd.to_numeric(freqs['A2_frequency'], errors='coerce')
candidates['Target_AF_A2'] = pd.to_numeric(candidates['Target_AF_A2'], errors='coerce')
candidates['Comp_AF_A2'] = pd.to_numeric(candidates['Comp_AF_A2'], errors='coerce')

# 1. Match status pie chart
plt.figure(figsize=(6, 6))
match_counts = match_status['matched_in_bim'].value_counts()
plt.pie(match_counts, labels=match_counts.index, autopct='%1.1f%%', colors=['#4CAF50', '#F44336'])
plt.title('Known Selection SNPs Match Status')
plt.savefig(os.path.join(FIG_DIR, 'pie_match_status_260522.png'), dpi=150, bbox_inches='tight')
plt.close()

# 2. Allele-frequency heatmap for known SNPs
plt.figure(figsize=(12, 8))
freqs_target = freqs[freqs['Population'].isin(TARGET_POPS)].copy()
pivot_freq = freqs_target.pivot(index="Population", columns="rsID", values="A2_frequency")
pivot_freq = pivot_freq.dropna(axis=1, how='all') # drop columns with all NaNs
sns.heatmap(pivot_freq, cmap="viridis", annot=True, fmt=".2f", cbar_kws={'label': 'A2 Frequency'})
plt.title('Known Selection SNPs Allele Frequency (Target Populations)')
plt.xticks(rotation=45, ha='right')
plt.savefig(os.path.join(FIG_DIR, 'heatmap_known_snps_260522.png'), dpi=150, bbox_inches='tight')
plt.close()

# 3. High-altitude adaptation SNP plot (EPAS1, EGLN1)
alt_snps = freqs_target[freqs_target['trait_or_adaptation'] == 'high-altitude adaptation'].copy()
if not alt_snps.empty:
    plt.figure(figsize=(10, 6))
    sns.barplot(data=alt_snps, x='Population', y='A2_frequency', hue='rsID')
    plt.title('High-Altitude Adaptation SNPs (EPAS1/EGLN1) in Target Populations')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Allele Frequency')
    plt.savefig(os.path.join(FIG_DIR, 'barplot_altitude_snps_260522.png'), dpi=150, bbox_inches='tight')
    plt.close()

# 4. Pigmentation-related SNP plot
pig_snps = freqs_target[freqs_target['trait_or_adaptation'].str.contains('pigmentation', na=False)].copy()
if not pig_snps.empty:
    plt.figure(figsize=(10, 6))
    sns.barplot(data=pig_snps, x='Population', y='A2_frequency', hue='rsID')
    plt.title('Pigmentation-related SNPs in Target Populations')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Allele Frequency')
    plt.savefig(os.path.join(FIG_DIR, 'barplot_pigmentation_snps_260522.png'), dpi=150, bbox_inches='tight')
    plt.close()

# 5. Missingness per known selection SNP
plt.figure(figsize=(10, 5))
miss_mean = freqs_target.groupby('rsID')['Missing_rate'].mean().reset_index()
sns.barplot(data=miss_mean, x='rsID', y='Missing_rate', color='salmon')
plt.title('Average Missingness Rate per Known Selection SNP (Targets)')
plt.xticks(rotation=45, ha='right')
plt.ylabel('Average Missingness')
plt.savefig(os.path.join(FIG_DIR, 'barplot_missingness_260522.png'), dpi=150, bbox_inches='tight')
plt.close()

# 6. Sample-size summary
# Total_alleles / 2 gives N called
freqs_target['N_called'] = freqs_target['Total_alleles'] / 2
n_mean = freqs_target.groupby('Population')['N_called'].max().reset_index()
plt.figure(figsize=(10, 5))
sns.barplot(data=n_mean, x='Population', y='N_called', color='skyblue')
plt.title('Max Called Sample Size per Target Population')
plt.xticks(rotation=45, ha='right')
plt.ylabel('Number of Individuals')
plt.savefig(os.path.join(FIG_DIR, 'barplot_sample_size_260522.png'), dpi=150, bbox_inches='tight')
plt.close()

# 7. Scatter plot of Target AF vs Comp AF for candidates
plt.figure(figsize=(8, 8))
sns.scatterplot(data=candidates, x='Comp_AF_A2', y='Target_AF_A2', hue='Category', s=60, alpha=0.8)
plt.plot([0, 1], [0, 1], 'k--', zorder=0)
plt.title('Candidate SNPs: Target vs Comparison Allele Frequencies')
plt.xlabel('Comparison AF')
plt.ylabel('Target AF')
plt.savefig(os.path.join(FIG_DIR, 'scatter_target_vs_comp_260522.png'), dpi=150, bbox_inches='tight')
plt.close()

# 8. Top 10 differentiated candidate SNPs plot
candidates['abs_diff'] = (candidates['Target_AF_A2'] - candidates['Comp_AF_A2']).abs()
top_diff = candidates.sort_values(by='abs_diff', ascending=False).head(10)
top_diff_melt = pd.melt(top_diff, id_vars=['rsID'], value_vars=['Target_AF_A2', 'Comp_AF_A2'], var_name='Group', value_name='AF')
plt.figure(figsize=(10, 6))
sns.barplot(data=top_diff_melt, x='rsID', y='AF', hue='Group')
plt.title('Top 10 Most Differentiated Candidate SNPs')
plt.xticks(rotation=45, ha='right')
plt.ylabel('Allele Frequency')
plt.savefig(os.path.join(FIG_DIR, 'barplot_top10_differentiated_260522.png'), dpi=150, bbox_inches='tight')
plt.close()

# 9. Fixed candidate counts per target pop
fixed = candidates[candidates['Fixed_Target_Pops'].notna() & (candidates['Fixed_Target_Pops'] != "")]
if not fixed.empty:
    fixed_counts = collections.Counter()
    for pops in fixed['Fixed_Target_Pops']:
        for p in pops.split(';'):
            fixed_counts[p] += 1
    fixed_df = pd.DataFrame(fixed_counts.items(), columns=['Population', 'Fixed_Count'])
    plt.figure(figsize=(10, 5))
    sns.barplot(data=fixed_df, x='Population', y='Fixed_Count', color='purple')
    plt.title('Number of Fixed Candidate SNPs per Target Population')
    plt.xticks(rotation=45, ha='right')
    plt.savefig(os.path.join(FIG_DIR, 'barplot_fixed_counts_260522.png'), dpi=150, bbox_inches='tight')
    plt.close()
else:
    # create empty placeholder
    fig, ax = plt.subplots(figsize=(6,4))
    ax.text(0.5, 0.5, "No fixed candidates found", ha='center', va='center')
    plt.savefig(os.path.join(FIG_DIR, 'barplot_fixed_counts_260522.png'), dpi=150)
    plt.close()

# 10. Characteristic SNP Summary Image
# This image combines the EPAS1 rs13419896 frequency and the top differentiated candidate.
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle('Characteristic SNP Findings Summary', fontsize=16, fontweight='bold')

# Panel 1: Known EPAS1 SNP if matched
epas1 = alt_snps[alt_snps['rsID'] == 'rs13419896']
if not epas1.empty:
    sns.barplot(data=epas1, x='Population', y='A2_frequency', ax=axes[0], color='dodgerblue')
    axes[0].set_title('High-Altitude Adaptation: EPAS1 (rs13419896)')
    axes[0].set_ylabel('Allele Frequency')
    axes[0].tick_params(axis='x', rotation=45)
else:
    axes[0].text(0.5, 0.5, "EPAS1 rs13419896 not matched", ha='center', va='center')

# Panel 2: Top Differentiated Candidate
if not top_diff.empty:
    top_snp = top_diff.iloc[0]
    axes[1].bar(['Target', 'Comparison'], [top_snp['Target_AF_A2'], top_snp['Comp_AF_A2']], color=['crimson', 'gray'])
    axes[1].set_title(f"Top Differentiated Novel Candidate: {top_snp['rsID']}")
    axes[1].set_ylabel('Allele Frequency')
    axes[1].text(0.5, 0.9, f"Target AF: {top_snp['Target_AF_A2']:.2f}\nComp AF: {top_snp['Comp_AF_A2']:.2f}", ha='center', va='center', transform=axes[1].transAxes, bbox=dict(facecolor='white', alpha=0.8))
else:
    axes[1].text(0.5, 0.5, "No novel candidates found", ha='center', va='center')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'characteristic_snp_summary_260522.png'), dpi=200, bbox_inches='tight')
plt.close()

print("Successfully generated 10 visualization files in the figures directory.")
