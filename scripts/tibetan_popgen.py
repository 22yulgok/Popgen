import os
import re
import sys
import warnings
import numpy as np
import pandas as pd
import scipy.linalg as la
import scipy.cluster.hierarchy as sch
from scipy.stats import chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

warnings.filterwarnings('ignore')

# ── Style & Config ─────────────────────────────────────────────────────────────
plt.rcParams.update({
    'figure.dpi': 150,
    'font.family': 'DejaVu Sans',
    'font.size': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.grid': True,
    'grid.alpha': 0.15,
    'grid.linestyle': '--',
    'axes.titlesize': 11,
    'axes.titleweight': 'bold',
    'axes.labelsize': 9,
    'legend.fontsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
})

# ── Paths ──────────────────────────────────────────────────────────────────────
PREFIX = '/home/duri_bae/popgen_lab/raw/Tibetan.fstat.1240K.260508'
FAM_PATH = PREFIX + '.fam'
BIM_PATH = PREFIX + '.bim'
BED_PATH = PREFIX + '.bed'
TARGET_META_PATH = '/home/duri_bae/popgen_lab/output/target_metadata_260522.csv'
SAMPLE_INFO_PATH = '/home/duri_bae/popgen_lab/output/SampleInfo_new_260522.csv'
AADR_PATH = '/home/duri_bae/popgen_lab/raw/v66.1240K.aadr.PUB.csv'
OUT_DIR = '/home/duri_bae/popgen_lab/output/rawvisual/Tib'
DATE_STAMP = '260522'

os.makedirs(OUT_DIR, exist_ok=True)

def save_fig(name):
    p = os.path.join(OUT_DIR, f'{name}_{DATE_STAMP}.png')
    plt.savefig(p, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  → Saved figure: {p}")
    return p

# ── Target Pop List & Color Palette ────────────────────────────────────────────
TARGET_POPS = [
    'GBSL', 'GBSL_old1', 'GBSL_old2_1', 'GBSL_old2_2', 'GBSL_old2_3',
    'Zongri5.1k', 'Zongri4.7k', 'Zongri4.5k', 'Zongri4.1k',
    'Yushu2.8k', 'Shannan3k', 'Chamdo2.8k_1',
    'Mbc4.4k_G1', 'Mbc4.4k_G2', 'Mbc4.4k_G0', 'Mbc4k_G2', 'Mbc3.5k_G0', 'Mbc3.5k_outlier',
    'SUI001.WGS'
]

# Note: Added Xingyi_EN here to avoid KeyError
REF_POPS = [
    'Mbuti.DG', 'Ami.DG', 'Atayal.DG', 'Haneyi', 'Boisman_MN', 'Baikal_EN',
    'DevilsCave_N', 'Upper_YR_LN', 'Yangshao_UYR', 'YR_LN', 'YR_MN', 'Chokhopani',
    'Anatolia_N', 'CHG', 'EHG', 'WHG', 'Yamnaya_Samara', 'Yana_UP', 'Onge.DG',
    'Xingyi_EN'
]

TARGET_COLORS = {
    'GBSL': '#2b5c8f', 'GBSL_old1': '#3a7cb8', 'GBSL_old2_1': '#5b9cd8', 'GBSL_old2_2': '#7bbcf8', 'GBSL_old2_3': '#9bdcf8',
    'Zongri5.1k': '#8b2b2b', 'Zongri4.7k': '#b83a3a', 'Zongri4.5k': '#d85b5b', 'Zongri4.1k': '#f87b7b',
    'Yushu2.8k': '#d4a373', 'Shannan3k': '#e9d8a6', 'Chamdo2.8k_1': '#ee9b00',
    'Mbc4.4k_G1': '#005f73', 'Mbc4.4k_G2': '#0a9396', 'Mbc4.4k_G0': '#94d2bd',
    'Mbc4k_G2': '#e07a5f', 'Mbc3.5k_G0': '#3d405b', 'Mbc3.5k_outlier': '#f4f1de',
    'SUI001.WGS': '#8338ec'
}

REF_COLORS = {
    'Mbuti.DG': '#7f7f7f', 'Ami.DG': '#1f77b4', 'Atayal.DG': '#aec7e8', 'Haneyi': '#ff7f0e',
    'Boisman_MN': '#2ca02c', 'Baikal_EN': '#98df8a', 'DevilsCave_N': '#d62728',
    'Upper_YR_LN': '#ff9896', 'Yangshao_UYR': '#9467bd', 'YR_LN': '#c5b0d5', 'YR_MN': '#8c564b',
    'Chokhopani': '#bcbd22', 'Anatolia_N': '#e377c2', 'CHG': '#f7b6d2', 'EHG': '#17becf',
    'WHG': '#9edae5', 'Yamnaya_Samara': '#dbdb8d', 'Yana_UP': '#c7c7c7', 'Onge.DG': '#db5f57',
    'Xingyi_EN': '#d62728'
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. Genotype Parsing & Loading
# ══════════════════════════════════════════════════════════════════════════════
print("Loading FAM and BIM files...")
fam = pd.read_csv(FAM_PATH, sep=r'\s+', header=None,
                  names=['FID','IID','PID','MID','SEX','PHENO'])
bim = pd.read_csv(BIM_PATH, sep=r'\t', header=None,
                  names=['CHR','SNP','CM','POS','A1','A2'])

N_SAM = len(fam)
N_SNP = len(bim)

# Map GBSL_old2 target sub-individuals explicitly
for idx, row in fam.iterrows():
    iid = row['IID']
    if iid == 'GBSL13':
        fam.at[idx, 'FID'] = 'GBSL_old2_1'
    elif iid == 'GBSL14':
        fam.at[idx, 'FID'] = 'GBSL_old2_2'
    elif iid == 'GBSL19':
        fam.at[idx, 'FID'] = 'GBSL_old2_3'

# Filter to keep only target and reference samples
selected_samples = fam[fam['FID'].isin(TARGET_POPS + REF_POPS)].copy()
selected_indices = selected_samples.index.tolist()
N_SEL = len(selected_samples)
print(f"Selected {N_SEL} samples out of {N_SAM} total (Targets: {selected_samples['FID'].isin(TARGET_POPS).sum()}, References: {selected_samples['FID'].isin(REF_POPS).sum()})")

# Load Target Metadata for details and dating
df_target = pd.read_csv(TARGET_META_PATH)
df_sample_info = pd.read_csv(SAMPLE_INFO_PATH)

# Create unique lookup for targets for geographic/temporal details
target_meta_map = {}
for _, row in df_target.iterrows():
    iid = str(row['IID']).strip()
    fid = str(row['FID']).strip()
    target_meta_map[iid] = row.to_dict()
    # Map GBSL_old2 sub-individuals
    if iid == 'GBSL13': target_meta_map[iid]['FID'] = 'GBSL_old2_1'
    elif iid == 'GBSL14': target_meta_map[iid]['FID'] = 'GBSL_old2_2'
    elif iid == 'GBSL19': target_meta_map[iid]['FID'] = 'GBSL_old2_3'

# Build a map of sample ages
aadr = pd.read_csv(AADR_PATH, dtype=str, low_memory=False)
aadr_date_numeric = pd.to_numeric(aadr[aadr.columns[10]], errors='coerce')
aadr_iid_map = dict(zip(aadr[aadr.columns[0]].str.strip(), aadr_date_numeric))
aadr_iid_map_2 = dict(zip(aadr[aadr.columns[2]].str.strip(), aadr_date_numeric))

def get_age_bp(row):
    iid = str(row['IID']).strip()
    fid = str(row['FID']).strip()
    pid = str(row['PID']).strip()
    
    # Check manual override for targets first
    if iid in target_meta_map:
        cal = str(target_meta_map[iid].get('Cal', ''))
        # Try PID suffix years (Mbc samples)
        m = re.search(r'(\d+\.?\d*)k', pid, re.IGNORECASE)
        if m:
            return float(m.group(1)) * 1000.0
        m = re.search(r'(\d+\.?\d*)k', fid, re.IGNORECASE)
        if m:
            return float(m.group(1)) * 1000.0
        # Check standard Cal text regex for BP
        m = re.search(r'(\d+)\s*BP', cal)
        if m:
            return float(m.group(1))
        m = re.search(r'mean:\s*(\d+)', cal)
        if m:
            return float(m.group(1))
        m = re.search(r'(\d+)-(\d+)\s*BCE', cal)
        if m:
            return (float(m.group(1)) + float(m.group(2))) / 2.0 + 1950.0
        
    # Check AADR mapping
    if iid in aadr_iid_map and not pd.isna(aadr_iid_map[iid]):
        return aadr_iid_map[iid]
    if iid in aadr_iid_map_2 and not pd.isna(aadr_iid_map_2[iid]):
        return aadr_iid_map_2[iid]
        
    # Fallback default values
    if 'Zongri5.1k' in fid: return 5100.0
    if 'Zongri4.7k' in fid: return 4700.0
    if 'Zongri4.5k' in fid: return 4500.0
    if 'Zongri4.1k' in fid: return 4100.0
    if 'Yushu2.8k' in fid: return 2800.0
    if 'Shannan3k' in fid: return 3000.0
    if 'Chamdo2.8k' in fid: return 2800.0
    if 'GBSL_old2' in fid: return 2550.0
    if 'GBSL_old1' in fid: return 2200.0
    if 'GBSL' in fid: return 1800.0
    if 'SUI001.WGS' in fid: return 3375.0
    
    # Try parsing text k suffix anywhere
    for txt in [fid, pid, iid]:
        m = re.search(r'(\d+\.?\d*)k', txt, re.IGNORECASE)
        if m:
            return float(m.group(1)) * 1000.0
            
    return np.nan

selected_samples['age_BP'] = selected_samples.apply(get_age_bp, axis=1)

# Memory-efficient BED parsing
CHUNK_SIZE = 10000
BYTES_PER_SNP = -(-N_SAM // 4)

lut = np.zeros((256, 4), dtype=np.uint8)
for byte in range(256):
    for bit_pair in range(4):
        lut[byte, bit_pair] = (byte >> (2 * bit_pair)) & 3

# PLINK BED code mapping to dosage (00->0, 01->NaN, 10->1, 11->2)
dosage_lut = np.array([0.0, np.nan, 1.0, 2.0], dtype=np.float32)

print("Decoding PLINK BED file for selected samples...")
G_list = []
with open(BED_PATH, 'rb') as f:
    magic = f.read(3)
    assert magic == b'\x6c\x1b\x01', f"Unexpected magic code: {magic.hex()}"
    
    snp_idx = 0
    while snp_idx < N_SNP:
        chunk_end = min(snp_idx + CHUNK_SIZE, N_SNP)
        chunk_n = chunk_end - snp_idx
        
        raw = np.frombuffer(f.read(chunk_n * BYTES_PER_SNP), dtype=np.uint8)
        raw = raw.reshape(chunk_n, BYTES_PER_SNP)
        
        # Decode chunk
        decoded = lut[raw]
        decoded = decoded.reshape(chunk_n, -1)[:, :N_SAM]
        
        # Keep only selected samples and convert to dosage float32
        sliced = dosage_lut[decoded[:, selected_indices]]
        G_list.append(sliced)
        
        snp_idx = chunk_end
        if snp_idx % 100000 == 0 or snp_idx == N_SNP:
            print(f"  Parsed {snp_idx:,} / {N_SNP:,} SNPs")

# Genotype matrix of shape (N_SNPs, N_selected_samples)
G = np.vstack(G_list)
print(f"Genotype matrix shape: {G.shape}")

# Precompute Block Jackknife block assignments
# Assign SNPs to 5 Mb blocks
block_size = 5_000_000
bim['block_id'] = bim['CHR'].astype(np.int32) * 1000 + (bim['POS'] // block_size).astype(np.int32)
unique_blocks = bim['block_id'].unique()
bim['block_idx'] = bim['block_id'].map({b: idx for idx, b in enumerate(unique_blocks)})
N_BLOCKS = len(unique_blocks)
print(f"Genome divided into {N_BLOCKS} blocks for Block Jackknife")

# ══════════════════════════════════════════════════════════════════════════════
# Sample Missingness and Call Rate Plots
# ══════════════════════════════════════════════════════════════════════════════
selected_samples['missing'] = np.isnan(G).sum(axis=0)
selected_samples['callrate'] = 1.0 - selected_samples['missing'] / N_SNP

print("Generating Sample Count and Missingness figures [Fig 4, 5]...")
# Fig 4: sample_count_by_target
plt.figure(figsize=(8, 4))
target_counts = selected_samples[selected_samples['FID'].isin(TARGET_POPS)]['FID'].value_counts()
colors = [TARGET_COLORS.get(pop, '#888888') for pop in target_counts.index]
target_counts.plot(kind='bar', color=colors, edgecolor='none')
plt.title("Sample Counts Across Target Early Tibetan Plateau Groups")
plt.xlabel("Population Group")
plt.ylabel("Number of Samples")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
save_fig("sample_count_by_target")

# Fig 5: snp_missingness_by_target
plt.figure(figsize=(9, 4.5))
box_data = []
box_labels = []
box_colors = []
for pop in TARGET_POPS:
    pop_callrates = selected_samples[selected_samples['FID'] == pop]['callrate'].values
    if len(pop_callrates) > 0:
        box_data.append(pop_callrates)
        box_labels.append(pop)
        box_colors.append(TARGET_COLORS.get(pop, '#888888'))

bp = plt.boxplot(box_data, labels=box_labels, patch_artist=True)
for patch, color in zip(bp['boxes'], box_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
    patch.set_edgecolor('black')
plt.title("Genotype Call Rates Across Target Populations")
plt.xlabel("Target Population")
plt.ylabel("Genotype Call Rate")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
save_fig("snp_missingness_by_target")

# ══════════════════════════════════════════════════════════════════════════════
# PCA Analyses (Broad and Focused)
# ══════════════════════════════════════════════════════════════════════════════
print("Running SVD-based PCA projection...")

def run_pca_projection(ref_pop_list, broad_flag=True):
    # Get reference and target indices
    ref_mask = selected_samples['FID'].isin(ref_pop_list)
    target_mask = selected_samples['FID'].isin(TARGET_POPS)
    
    ref_idx_sub = np.where(ref_mask)[0]
    target_idx_sub = np.where(target_mask)[0]
    
    # Filter SNPs: missingness in references, MAF > 1%
    G_ref = G[:, ref_idx_sub]
    missing_rate_ref = np.isnan(G_ref).mean(axis=1)
    
    ref_mean = np.nanmean(G_ref, axis=1)
    ref_freq = ref_mean / 2.0
    ref_maf = np.minimum(ref_freq, 1.0 - ref_freq)
    
    # Loosen missingness threshold dynamically to get enough SNPs
    valid_snps = None
    for thresh in [0.2, 0.5, 0.7, 0.8, 0.9, 0.95]:
        valid_snps = (missing_rate_ref < thresh) & (ref_maf > 0.01)
        if valid_snps.sum() >= 5000:
            print(f"  Selected missingness threshold < {thresh} for PCA, yielding {valid_snps.sum():,} SNPs")
            break
            
    if valid_snps is None or valid_snps.sum() < 100:
        # Final fallback: use all SNPs with any variation
        valid_snps = ref_maf > 0.001
        print(f"  PCA Fallback: Selected all polymorphic SNPs ({valid_snps.sum():,})")
    
    G_ref_filtered = G_ref[valid_snps, :]
    ref_mean_filtered = ref_mean[valid_snps]
    ref_std_filtered = np.nanstd(G_ref_filtered, axis=1)
    
    # Remove standard deviations close to zero
    non_zero_std = ref_std_filtered > 1e-5
    G_ref_filtered = G_ref_filtered[non_zero_std, :]
    ref_mean_filtered = ref_mean_filtered[non_zero_std]
    ref_std_filtered = ref_std_filtered[non_zero_std]
    valid_snps_indices = np.where(valid_snps)[0][non_zero_std]
    
    # Standardize references and mean-impute remaining missingness
    G_ref_std = (G_ref_filtered - ref_mean_filtered[:, None]) / ref_std_filtered[:, None]
    G_ref_std[np.isnan(G_ref_std)] = 0.0
    
    # Run SVD on references
    U, S, Vt = la.svd(G_ref_std, full_matrices=False)
    
    # Reference coordinates
    ref_coords = G_ref_std.T @ U[:, :4]
    
    # Least-squares project targets to avoid missingness bias
    target_coords = []
    G_target = G[valid_snps_indices, :][:, target_idx_sub]
    
    for i in range(len(target_idx_sub)):
        g_sample = G_target[:, i]
        c = np.zeros(4)
        non_missing = ~np.isnan(g_sample)
        if non_missing.sum() >= 100:
            y = (g_sample[non_missing] - ref_mean_filtered[non_missing]) / ref_std_filtered[non_missing]
            X = U[non_missing, :4]
            c, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        target_coords.append(c)
        
    target_coords = np.array(target_coords)
    
    return ref_idx_sub, target_idx_sub, ref_coords, target_coords, valid_snps_indices, U[:, :4]

# Broad PCA
broad_refs = [pop for pop in REF_POPS if pop != 'Mbuti.DG']
ref_idx_b, target_idx_b, ref_coords_b, target_coords_b, snp_idx_b, U_b = run_pca_projection(broad_refs, broad_flag=True)

# Focused PCA
focused_refs = [pop for pop in REF_POPS if pop not in ['Anatolia_N', 'CHG', 'EHG', 'WHG', 'Yamnaya_Samara', 'Yana_UP', 'Mbuti.DG']]
ref_idx_f, target_idx_f, ref_coords_f, target_coords_f, snp_idx_f, U_f = run_pca_projection(focused_refs, broad_flag=False)

# Save PCA coordinates to CSV
pca_df = pd.DataFrame(columns=['IID', 'FID', 'Type', 'PC1', 'PC2', 'PC3', 'PC4'])
for i, idx in enumerate(ref_idx_b):
    row = selected_samples.iloc[idx]
    pca_df.loc[len(pca_df)] = [row['IID'], row['FID'], 'Reference', ref_coords_b[i, 0], ref_coords_b[i, 1], ref_coords_b[i, 2], ref_coords_b[i, 3]]
for i, idx in enumerate(target_idx_b):
    row = selected_samples.iloc[idx]
    pca_df.loc[len(pca_df)] = [row['IID'], row['FID'], 'Target', target_coords_b[i, 0], target_coords_b[i, 1], target_coords_b[i, 2], target_coords_b[i, 3]]
pca_df.to_csv(os.path.join(OUT_DIR, 'pca_loadings_260522.csv'), index=False)

# Plotting Broad PCA [Fig 1]
plt.figure(figsize=(8, 6))
# Plot references
for pop in REF_POPS:
    mask = selected_samples.iloc[ref_idx_b]['FID'] == pop
    if mask.sum() > 0:
        plt.scatter(ref_coords_b[mask, 0], ref_coords_b[mask, 1],
                    color=REF_COLORS.get(pop, '#888888'), label=pop, alpha=0.5, s=25, edgecolor='none')
# Plot targets prominently
for pop in TARGET_POPS:
    mask = selected_samples.iloc[target_idx_b]['FID'] == pop
    if mask.sum() > 0:
        plt.scatter(target_coords_b[mask, 0], target_coords_b[mask, 1],
                    color=TARGET_COLORS.get(pop, '#ff007f'), label=pop, s=60, edgecolor='black', zorder=10)
plt.title("Broad PCA (Projected Target Early Tibetan Populations)")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2)
plt.tight_layout()
save_fig("pca_broad_targets")

# Plotting Focused PCA [Fig 2]
plt.figure(figsize=(8, 6))
for pop in focused_refs:
    mask = selected_samples.iloc[ref_idx_f]['FID'] == pop
    if mask.sum() > 0:
        plt.scatter(ref_coords_f[mask, 0], ref_coords_f[mask, 1],
                    color=REF_COLORS.get(pop, '#888888'), label=pop, alpha=0.5, s=25, edgecolor='none')
for pop in TARGET_POPS:
    mask = selected_samples.iloc[target_idx_f]['FID'] == pop
    if mask.sum() > 0:
        plt.scatter(target_coords_f[mask, 0], target_coords_f[mask, 1],
                    color=TARGET_COLORS.get(pop, '#ff007f'), label=pop, s=60, edgecolor='black', zorder=10)
plt.title("Focused PCA on East Asian and Tibetan related Groups")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=1)
plt.tight_layout()
save_fig("pca_tibetan_focused")

# Zoomed-in Target Labels PCA [Fig 3]
plt.figure(figsize=(7, 7))
for pop in TARGET_POPS:
    mask = selected_samples.iloc[target_idx_f]['FID'] == pop
    if mask.sum() > 0:
        plt.scatter(target_coords_f[mask, 0], target_coords_f[mask, 1],
                    color=TARGET_COLORS.get(pop, '#ff007f'), s=80, edgecolor='black', zorder=10)
        # Add labels
        for idx_c, val in enumerate(mask):
            if val:
                row = selected_samples.iloc[target_idx_f[idx_c]]
                plt.text(target_coords_f[idx_c, 0] + 0.05, target_coords_f[idx_c, 1] + 0.05,
                         row['IID'], fontsize=7, alpha=0.8, weight='bold')
plt.title("Zoomed-in Target PCA with Individual Sample Labels")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
save_fig("pca_target_labels")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Block Jackknife Outgroup-f3 & f4 Computations
# ══════════════════════════════════════════════════════════════════════════════
print("Computing Outgroup-f3 and f4 statistics...")

# Precompute population frequencies and block-level statistics
pop_indices = {}
for pop in TARGET_POPS + REF_POPS:
    pop_indices[pop] = np.where(selected_samples['FID'] == pop)[0]

# Calculate population frequencies at all SNPs
pop_freqs = {}
for pop, idxs in pop_indices.items():
    if len(idxs) > 0:
        # Sum allele dosage / (2 * non_missing count)
        pop_freqs[pop] = np.nanmean(G[:, idxs], axis=1) / 2.0
    else:
        pop_freqs[pop] = np.full(N_SNP, np.nan)

# Jackknife blocks precomputation
def run_block_jackknife(num_vals, den_vals, bim):
    block_idx = bim['block_idx'].values
    
    # Remove nan-containing SNPs
    valid = ~np.isnan(num_vals) & ~np.isnan(den_vals)
    num_vals = num_vals[valid]
    den_vals = den_vals[valid]
    block_idx = block_idx[valid]
    
    if len(num_vals) < 10:
        return np.nan, np.nan, np.nan, 0
        
    total_num = np.sum(num_vals)
    total_den = np.sum(den_vals)
    overall_stat = total_num / total_den
    
    # Calculate sums for each block
    block_num_sums = np.bincount(block_idx, weights=num_vals, minlength=N_BLOCKS)
    block_den_sums = np.bincount(block_idx, weights=den_vals, minlength=N_BLOCKS)
    
    # Jackknife replicates
    jack_num = total_num - block_num_sums
    jack_den = total_den - block_den_sums
    
    valid_reps = jack_den != 0
    reps = jack_num[valid_reps] / jack_den[valid_reps]
    N = len(reps)
    
    if N < 2:
        return overall_stat, np.nan, np.nan, len(num_vals)
        
    # Jackknife standard error
    se = np.sqrt(((N - 1) / N) * np.sum((reps - reps.mean()) ** 2))
    z = overall_stat / se if se > 0 else np.nan
    
    return overall_stat, se, z, len(num_vals)

# Compute f3(Mbuti.DG; Target, Reference)
f3_results = []
outgroup = 'Mbuti.DG'
p_O = pop_freqs[outgroup]

print("Calculating Outgroup-f3 statistics...")
for target in TARGET_POPS:
    p_A = pop_freqs[target]
    if np.isnan(p_A).all():
        continue
    for ref in REF_POPS:
        if ref == outgroup:
            continue
        p_B = pop_freqs[ref]
        
        # Unbiased outgroup-f3: (p_A - p_O)*(p_B - p_O)
        num = (p_A - p_O) * (p_B - p_O)
        den = np.ones(N_SNP)
        
        val, se, z, n_snps = run_block_jackknife(num, den, bim)
        f3_results.append({
            'Target': target,
            'Reference': ref,
            'f3': val,
            'se': se,
            'z': z,
            'nSNPs': n_snps
        })

df_f3 = pd.DataFrame(f3_results)
df_f3.to_csv(os.path.join(OUT_DIR, 'f3_outgroup_summary_260522.csv'), index=False)

# Ranked Outgroup-f3 plots [Fig 6, 7, 8]
def plot_ranked_f3(target_list, name_str, filename):
    plt.figure(figsize=(8, 5))
    target_f3 = df_f3[df_f3['Target'].isin(target_list)].groupby('Reference')['f3'].mean().sort_values(ascending=False)
    target_f3_se = df_f3[df_f3['Target'].isin(target_list)].groupby('Reference')['se'].mean().loc[target_f3.index]
    
    colors = [REF_COLORS.get(ref, '#888888') for ref in target_f3.index]
    plt.barh(target_f3.index, target_f3.values, xerr=target_f3_se.values, color=colors, edgecolor='none', alpha=0.8)
    plt.title(f"Ranked Outgroup-f3 Affinity for {name_str}\n" + r"$f_3(\mathrm{Mbuti.DG}; \mathrm{Target}, \mathrm{Reference})$")
    plt.xlabel("f3 value")
    plt.tight_layout()
    save_fig(filename)

plot_ranked_f3(['Zongri5.1k', 'Zongri4.7k', 'Zongri4.5k', 'Zongri4.1k'], "Zongri Plateau Groups", "outgroup_f3_ranked_zongri")
plot_ranked_f3(['GBSL', 'GBSL_old1', 'GBSL_old2_1', 'GBSL_old2_2', 'GBSL_old2_3'], "Gebusailu (GBSL) Groups", "outgroup_f3_ranked_gbsl")
plot_ranked_f3(['Mbc4.4k_G1', 'Mbc4.4k_G2', 'Mbc4.4k_G0', 'Mbc4k_G2', 'Mbc3.5k_G0', 'Mbc3.5k_outlier'], "Mabu Co (Mbc) Groups", "outgroup_f3_ranked_mabuco")

# Fig 9: Heatmap of outgroup f3
f3_pivot = df_f3.pivot(index='Target', columns='Reference', values='f3').loc[TARGET_POPS]
plt.figure(figsize=(10, 6))
sns.heatmap(f3_pivot, cmap='coolwarm', annot=True, fmt=".4f", annot_kws={"size": 6})
plt.title("Outgroup-f3 Heatmap: Target Plateau Groups vs. References")
plt.tight_layout()
save_fig("f3_target_heatmap")

# Compute f4 Symmetry Tests
f4_results = []
print("Calculating f4 symmetry tests...")

# Setup f4 test lists: f4(Mbuti.DG, Target; Pop3, Pop4)
f4_tests = []
for target in TARGET_POPS:
    # Yellow river vs Southern/Southeast Asian ancestry
    f4_tests.append((target, 'Upper_YR_LN', 'Xingyi_EN'))
    # Internal comparison: Zongri vs Himalayan Chokhopani
    f4_tests.append((target, 'Zongri5.1k', 'Chokhopani'))

for target, p3, p4 in f4_tests:
    p_A = pop_freqs[target]
    p_B = pop_freqs[p3]
    p_C = pop_freqs[p4]
    
    if np.isnan(p_A).all() or np.isnan(p_B).all() or np.isnan(p_C).all():
        continue
        
    num = (p_O - p_A) * (p_B - p_C)
    den = np.ones(N_SNP)
    
    val, se, z, n_snps = run_block_jackknife(num, den, bim)
    f4_results.append({
        'Target': target,
        'Pop3': p3,
        'Pop4': p4,
        'f4': val,
        'se': se,
        'z': z,
        'nSNPs': n_snps
    })

df_f4 = pd.DataFrame(f4_results)
df_f4.to_csv(os.path.join(OUT_DIR, 'f4_symmetry_summary_260522.csv'), index=False)

# Forest plot helper
def plot_f4_forest(p3, p4, title, filename):
    sub_df = df_f4[(df_f4['Pop3'] == p3) & (df_f4['Pop4'] == p4)]
    if len(sub_df) == 0:
        return
    sub_df = sub_df.sort_values('f4')
    
    plt.figure(figsize=(7, 6))
    colors = ['#d73027' if z_val >= 3.0 else ('#4575b4' if z_val <= -3.0 else '#4d4d4d') for z_val in sub_df['z']]
    
    plt.errorbar(sub_df['f4'], sub_df['Target'], xerr=3 * sub_df['se'], fmt='o', color='lightgrey', elinewidth=1, capsize=2, label='3 SE')
    plt.errorbar(sub_df['f4'], sub_df['Target'], xerr=1.96 * sub_df['se'], fmt='o', color='grey', elinewidth=1.5, capsize=3, label='1.96 SE')
    plt.scatter(sub_df['f4'], sub_df['Target'], color=colors, s=30, zorder=10, label='f4 estimate')
    
    plt.axvline(0, color='black', linestyle='--', alpha=0.5)
    plt.title(title + f"\nf4(Mbuti.DG, Target; {p3}, {p4})")
    plt.xlabel("f4 value")
    plt.ylabel("Target Population")
    plt.legend()
    plt.tight_layout()
    save_fig(filename)

# Fig 10 & 11 Forest plots
plot_f4_forest('Upper_YR_LN', 'Xingyi_EN', "Ancestry Test: Yellow River vs. Xingyi EN", "f4_xingyi_vs_yellowriver")
plot_f4_forest('Zongri5.1k', 'Chokhopani', "Internal Structure: Zongri vs. Himalayan Chokhopani", "f4_internal_tibetan_comparison")

# Fig 12: f4 Z-score significance plot
plt.figure(figsize=(8, 5))
plt.scatter(df_f4['z'], df_f4['Target'] + " (" + df_f4['Pop3'] + ")", c=np.abs(df_f4['z']), cmap='OrRd', edgecolor='black', s=45, zorder=5)
plt.axvline(-3, color='#e63946', linestyle='--', alpha=0.8, label='|Z|=3 threshold')
plt.axvline(3, color='#e63946', linestyle='--')
plt.axvline(0, color='grey', linestyle='-', alpha=0.5)
plt.title("f4 Test Z-score Significance Summary")
plt.xlabel("Z-score")
plt.ylabel("Test Target & Contrast")
plt.legend()
plt.tight_layout()
save_fig("f4_zscore_significance")


# ══════════════════════════════════════════════════════════════════════════════
# 3. Native qpWave and qpAdm Implementation
# ══════════════════════════════════════════════════════════════════════════════
print("Running qpWave and qpAdm modeling...")

def get_f4_replicates(target, left_pops, right_pops):
    base_L = left_pops[-1]
    base_R = right_pops[0]
    
    L_test = [target] + left_pops[:-1]
    n_L = len(L_test)
    n_R = len(right_pops) - 1
    dim = n_L * n_R
    
    block_idx = bim['block_idx'].values
    
    block_nums = np.zeros((N_BLOCKS, dim))
    block_dens = np.zeros((N_BLOCKS, dim))
    
    p_baseL = pop_freqs[base_L]
    p_baseR = pop_freqs[base_R]
    
    idx_col = 0
    for l in L_test:
        p_L = pop_freqs[l]
        for r_idx in range(1, len(right_pops)):
            r = right_pops[r_idx]
            p_R = pop_freqs[r]
            
            num_val = (p_baseL - p_L) * (p_baseR - p_R)
            valid = ~np.isnan(num_val)
            
            num_val_clean = np.zeros(N_SNP)
            num_val_clean[valid] = num_val[valid]
            
            block_nums[:, idx_col] = np.bincount(block_idx, weights=num_val_clean, minlength=N_BLOCKS)
            block_dens[:, idx_col] = np.bincount(block_idx, weights=valid.astype(np.float64), minlength=N_BLOCKS)
            idx_col += 1
            
    total_num = block_nums.sum(axis=0)
    total_den = block_dens.sum(axis=0)
    
    # Avoid division by zero
    total_den_safe = np.where(total_den == 0, 1.0, total_den)
    overall_stat = total_num / total_den_safe
    
    reps = []
    for b in range(N_BLOCKS):
        jack_num = total_num - block_nums[b]
        jack_den = total_den - block_dens[b]
        jack_den_safe = np.where(jack_den == 0, 1.0, jack_den)
        reps.append(jack_num / jack_den_safe)
    reps = np.array(reps)
    
    mean_reps = reps.mean(axis=0)
    V = ((N_BLOCKS - 1) / N_BLOCKS) * ((reps - mean_reps).T @ (reps - mean_reps))
    
    return overall_stat, V, reps, int(total_den.mean())

def run_qpadm(target, sources, right_pops):
    k = len(sources)
    o = len(right_pops)
    
    stat, V_z, reps, n_snps = get_f4_replicates(target, sources, right_pops)
    
    y = stat[:o-1]
    A = stat[o-1:].reshape(k-1, o-1).T
    
    alpha = np.ones(k-1) / (k-1)
    
    # Regularize V_z to avoid singularity issues in matrix inversion
    V_z_reg = V_z + np.eye(V_z.shape[0]) * 1e-7
    
    for itr in range(10):
        C = np.zeros((o-1, (o-1)*k))
        for j in range(o-1):
            C[j, j] = 1.0
            for a in range(k-1):
                C[j, (a+1)*(o-1) + j] = -alpha[a]
                
        Sigma = C @ V_z_reg @ C.T
        # Add small ridge to Sigma to ensure stability
        Sigma += np.eye(Sigma.shape[0]) * 1e-8
        Sigma_inv = la.pinv(Sigma)
        
        # GLS Estimate
        cov_alpha = la.pinv(A.T @ Sigma_inv @ A)
        alpha = cov_alpha @ A.T @ Sigma_inv @ y
        
    alpha_base = 1.0 - np.sum(alpha)
    proportions = np.append(alpha, alpha_base)
    
    se_alpha = np.sqrt(np.maximum(np.diag(cov_alpha), 0.0))
    se_base = np.sqrt(np.maximum(np.sum(cov_alpha), 0.0))
    standard_errors = np.append(se_alpha, se_base)
    
    chi2_val = (y - A @ alpha).T @ Sigma_inv @ (y - A @ alpha)
    df = o - k
    p_val = 1.0 - chi2.cdf(chi2_val, df) if df > 0 else 1.0
    
    return proportions, standard_errors, p_val, chi2_val, df, n_snps

# qpWave Rank 0 Clade Heatmap
print("Running qpWave Clade tests...")
qpwave_matrix = np.zeros((len(TARGET_POPS), len(TARGET_POPS)))
for i, target_i in enumerate(TARGET_POPS):
    for j, target_j in enumerate(TARGET_POPS):
        if i == j:
            qpwave_matrix[i, j] = 1.0
            continue
        try:
            stat, V, _, _ = get_f4_replicates(target_i, [target_j], REF_POPS)
            # Regularize V
            V_reg = V + np.eye(V.shape[0]) * 1e-7
            chi2_val = stat.T @ la.pinv(V_reg) @ stat
            df = len(REF_POPS) - 1
            p_val = 1.0 - chi2.cdf(chi2_val, df)
            qpwave_matrix[i, j] = p_val
        except Exception:
            qpwave_matrix[i, j] = np.nan

# Plotting qpWave Rank Heatmap [Fig 13]
plt.figure(figsize=(10, 8))
sns.heatmap(qpwave_matrix, xticklabels=TARGET_POPS, yticklabels=TARGET_POPS,
            cmap='viridis', annot=True, fmt=".2f", annot_kws={"size": 5})
plt.title("qpWave Clade Test Rank 0 p-values\n(High p-value indicates consistency with being a clade)")
plt.tight_layout()
save_fig("qpwave_rank_summary")

# Run qpAdm models
print("Running qpAdm admixture modeling...")
qpadm_results = []
models = [
    # GBSL targets
    ('GBSL', ['Chokhopani', 'Upper_YR_LN']),
    ('GBSL_old1', ['Chokhopani', 'Upper_YR_LN']),
    ('GBSL_old2_1', ['Chokhopani', 'Upper_YR_LN']),
    ('GBSL_old2_2', ['Chokhopani', 'Upper_YR_LN']),
    ('GBSL_old2_3', ['Chokhopani', 'Upper_YR_LN']),
    # Zongri targets
    ('Zongri4.5k', ['Zongri5.1k', 'Upper_YR_LN']),
    ('Zongri4.1k', ['Zongri5.1k', 'Upper_YR_LN']),
    # Mbc targets
    ('Mbc4.4k_G2', ['Mbc3.5k_outlier', 'Upper_YR_LN']),
    ('Mbc4k_G2', ['Mbc3.5k_outlier', 'Upper_YR_LN']),
    ('Mbc3.5k_G0', ['Mbc3.5k_outlier', 'Upper_YR_LN']),
    # SUI001.WGS target
    ('SUI001.WGS', ['Chokhopani', 'Upper_YR_LN'])
]

right_set = [ref for ref in REF_POPS if ref not in ['Upper_YR_LN', 'Chokhopani', 'Mbc3.5k_outlier', 'Zongri5.1k']]

for target, sources in models:
    try:
        props, ses, pval, chi_stat, df, n_snps = run_qpadm(target, sources, right_set)
        qpadm_results.append({
            'Target': target,
            'Source1': sources[0],
            'Source2': sources[1],
            'Prop1': props[0],
            'Prop2': props[1],
            'SE1': ses[0],
            'SE2': ses[1],
            'p_value': pval,
            'Chi2': chi_stat,
            'nSNPs': n_snps
        })
    except Exception as e:
        print(f"Failed to fit model {target} -> {sources}: {e}")

df_qpadm = pd.DataFrame(qpadm_results)
df_qpadm.to_csv(os.path.join(OUT_DIR, 'qpadm_model_summary_260522.csv'), index=False)

# qpAdm model pvalues bar plot [Fig 14]
plt.figure(figsize=(8, 4))
plt.bar(df_qpadm['Target'], df_qpadm['p_value'], color='#2a9d8f', edgecolor='black', alpha=0.8)
plt.axhline(0.05, color='red', linestyle='--', label='p=0.05 threshold')
plt.title("qpAdm Model Fit p-values (2-source models)")
plt.xlabel("Target Population")
plt.ylabel("Fit p-value")
plt.xticks(rotation=45, ha='right')
plt.legend()
plt.tight_layout()
save_fig("qpadm_model_pvalues")

# qpAdm admixture proportions stacked bar plot [Fig 15]
plt.figure(figsize=(9, 5))
ind = np.arange(len(df_qpadm))
width = 0.5

p1 = plt.bar(ind, df_qpadm['Prop1'], width, color='#3a86c8', edgecolor='none', label='Source 1 (Plateau/Himalaya)')
p2 = plt.bar(ind, df_qpadm['Prop2'], width, bottom=df_qpadm['Prop1'], color='#ff9f1c', edgecolor='none', label='Source 2 (Yellow River)')

plt.errorbar(ind, df_qpadm['Prop1'], yerr=df_qpadm['SE1'], fmt='none', ecolor='black', elinewidth=1.5, capsize=3)

plt.title("qpAdm Admixture Proportions with ±1 SE")
plt.xticks(ind, df_qpadm['Target'], rotation=45, ha='right')
plt.ylabel("Ancestry Proportions")
plt.ylim(0, 1.1)
plt.legend()
plt.tight_layout()
save_fig("qpadm_admixture_proportions")


# ══════════════════════════════════════════════════════════════════════════════
# 4. Hudson Fst Pairwise Matrix & Clustering Trees
# ══════════════════════════════════════════════════════════════════════════════
print("Calculating Hudson Fst matrix...")

def compute_hudson_fst(popA, popB):
    pA = pop_freqs[popA]
    pB = pop_freqs[popB]
    
    nA = 2 * len(pop_indices[popA])
    nB = 2 * len(pop_indices[popB])
    
    if nA < 2 or nB < 2:
        return np.nan
        
    num = (pA - pB)**2 - (pA*(1.0 - pA))/(nA - 1) - (pB*(1.0 - pB))/(nB - 1)
    den = pA*(1.0 - pB) + pB*(1.0 - pA)
    
    valid = ~np.isnan(num) & ~np.isnan(den) & (den != 0)
    if valid.sum() < 100:
        return np.nan
        
    return np.sum(num[valid]) / np.sum(den[valid])

# Calculate pairwise Fst for targets and references
fst_matrix = np.zeros((len(TARGET_POPS), len(TARGET_POPS)))
for i, pop1 in enumerate(TARGET_POPS):
    for j, pop2 in enumerate(TARGET_POPS):
        if i == j:
            fst_matrix[i, j] = 0.0
        elif i > j:
            val = compute_hudson_fst(pop1, pop2)
            fst_matrix[i, j] = val
            fst_matrix[j, i] = val

df_fst = pd.DataFrame(fst_matrix, index=TARGET_POPS, columns=TARGET_POPS)
df_fst.to_csv(os.path.join(OUT_DIR, 'fst_target_matrix_260522.csv'))

# Fig 16: fst_heatmap
plt.figure(figsize=(9, 7.5))
sns.heatmap(df_fst, cmap='Reds', annot=True, fmt=".4f", annot_kws={"size": 6})
plt.title("Pairwise Hudson Fst Matrix Between Target Populations")
plt.tight_layout()
save_fig("fst_heatmap")

# Fig 17: hierarchical clustering dendrogram based on Fst
plt.figure(figsize=(8, 5))
df_fst_clean = df_fst.fillna(0.2)
Z_linkage = sch.linkage(sch.distance.squareform(df_fst_clean), method='average')
sch.dendrogram(Z_linkage, labels=TARGET_POPS, orientation='left', color_threshold=0.1)
plt.title("Neighbor-Joining/UPGMA Hierarchical Tree (Hudson Fst Distance)")
plt.xlabel("Fst Genetic Distance")
plt.tight_layout()
save_fig("fst_nj_tree")


# ══════════════════════════════════════════════════════════════════════════════
# 5. Individual-level IBS Sharing, MDS & Heatmap
# ══════════════════════════════════════════════════════════════════════════════
print("Calculating individual IBS sharing matrix...")
targets_only = selected_samples[selected_samples['FID'].isin(TARGET_POPS)]
N_T = len(targets_only)
G_T = G[:, targets_only.index]

ibs_matrix = np.zeros((N_T, N_T))
for i in range(N_T):
    g_i = G_T[:, i]
    for j in range(i, N_T):
        g_j = G_T[:, j]
        diff = np.abs(g_i - g_j)
        valid = ~np.isnan(diff)
        if valid.sum() > 100:
            sharing = 1.0 - (np.sum(diff[valid]) / (2.0 * valid.sum()))
            ibs_matrix[i, j] = sharing
            ibs_matrix[j, i] = sharing
        else:
            ibs_matrix[i, j] = np.nan
            ibs_matrix[j, i] = np.nan

# Convert to distance matrix
ibs_dist = 1.0 - ibs_matrix

df_ibs = pd.DataFrame(ibs_matrix, index=targets_only['IID'], columns=targets_only['IID'])
df_ibs.to_csv(os.path.join(OUT_DIR, 'ibs_individual_matrix_260522.csv'))

# Fig 18: ibs_individual_heatmap
plt.figure(figsize=(10, 8.5))
sns.heatmap(df_ibs, cmap='YlGnBu_r', annot=False)
plt.title("Individual-level Identity by State (IBS) Sharing Matrix")
plt.tight_layout()
save_fig("ibs_individual_heatmap")

# Fig 19: classical MDS on IBS distance
print("Running Classical MDS Analysis...")
ibs_dist_clean = np.nan_to_num(ibs_dist, nan=np.nanmean(ibs_dist))
H = np.eye(N_T) - np.ones((N_T, N_T)) / N_T
B = -0.5 * (H @ (ibs_dist_clean ** 2) @ H)
evals, evecs = la.eigh(B)
idx = np.argsort(evals)[::-1]
evals = evals[idx]
evecs = evecs[:, idx]
mds_coords = evecs[:, :2] * np.sqrt(np.maximum(evals[:2], 0.0))

plt.figure(figsize=(8, 7))
for pop in TARGET_POPS:
    mask = targets_only['FID'] == pop
    if mask.sum() > 0:
        plt.scatter(mds_coords[mask, 0], mds_coords[mask, 1],
                    color=TARGET_COLORS.get(pop, '#ff007f'), label=pop, s=70, edgecolor='black', zorder=10)
        
        # Label sample IDs
        for idx_c, val in enumerate(mask):
            if val:
                plt.text(mds_coords[idx_c, 0] + 0.002, mds_coords[idx_c, 1] + 0.002,
                         targets_only.iloc[idx_c]['IID'], fontsize=6, alpha=0.8)
                
plt.title("Multidimensional Scaling (MDS) on IBS sharing distance")
plt.xlabel("MDS Dimension 1")
plt.ylabel("MDS Dimension 2")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
save_fig("ibs_mds")


# ══════════════════════════════════════════════════════════════════════════════
# 6. Temporal Cline & Geographical Timeline
# ══════════════════════════════════════════════════════════════════════════════
print("Modeling temporal genetic clines...")
fst_to_yr = []
for pop in TARGET_POPS:
    val = compute_hudson_fst(pop, 'Upper_YR_LN')
    fst_to_yr.append(val)
    
targets_only['Fst_to_YR'] = targets_only['FID'].map(dict(zip(TARGET_POPS, fst_to_yr)))

# Fig 20: temporal_cline (Distance to YR vs BP age)
plt.figure(figsize=(8, 5.5))
targets_with_age = targets_only.dropna(subset=['age_BP', 'Fst_to_YR'])

for pop in TARGET_POPS:
    mask = targets_with_age['FID'] == pop
    if mask.sum() > 0:
        plt.scatter(targets_with_age[mask]['age_BP'], targets_with_age[mask]['Fst_to_YR'],
                    color=TARGET_COLORS.get(pop, '#888888'), label=pop, s=80, edgecolor='black', zorder=10)

# Regression line
if len(targets_with_age) > 2:
    from scipy.stats import linregress
    res = linregress(targets_with_age['age_BP'], targets_with_age['Fst_to_YR'])
    x_range = np.linspace(targets_with_age['age_BP'].min() - 200, targets_with_age['age_BP'].max() + 200, 100)
    plt.plot(x_range, res.intercept + res.slope * x_range, color='grey', linestyle='--', alpha=0.8,
             label=f'R={res.rvalue:.3f}, p={res.pvalue:.4f}')

plt.title("Genetic Cline Over Time (Fst Distance to Upper_YR_LN vs BP age)")
plt.xlabel("Sample calibrated Age (BP)")
plt.ylabel("Hudson Fst Genetic Distance to Upper_YR_LN")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.gca().invert_xaxis()
plt.tight_layout()
save_fig("temporal_cline")


# Fig 21: Geographic-Temporal Timeline
print("Generating timeline visualization...")
targets_only['latitude'] = targets_only['IID'].map(lambda x: float(target_meta_map.get(x, {}).get('latitude', 30.0)) if x in target_meta_map and target_meta_map[x].get('latitude', '') != 'unknown' else 30.0)

plt.figure(figsize=(9, 6))
targets_timeline = targets_only.dropna(subset=['age_BP', 'latitude'])

for pop in TARGET_POPS:
    mask = targets_timeline['FID'] == pop
    if mask.sum() > 0:
        sizes = targets_timeline[mask]['callrate'] * 150
        plt.scatter(targets_timeline[mask]['age_BP'], targets_timeline[mask]['latitude'],
                    color=TARGET_COLORS.get(pop, '#888888'), label=pop, s=sizes, edgecolor='black', alpha=0.85, zorder=5)

plt.title("Geographic-Temporal Timeline of Target early Tibetan Plateau Groups")
plt.xlabel("Calibrated Age (BP)")
plt.ylabel("Geographic Latitude (°N)")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.gca().invert_xaxis()
plt.tight_layout()
save_fig("target_timeline")

# Save other summary CSV files as requested
df_f4.to_csv(os.path.join(OUT_DIR, 'f4_symmetry_summary_260522.csv'), index=False)
df_qpadm.to_csv(os.path.join(OUT_DIR, 'qpadm_model_summary_260522.csv'), index=False)
df_fst.to_csv(os.path.join(OUT_DIR, 'fst_target_matrix_260522.csv'))
df_ibs.to_csv(os.path.join(OUT_DIR, 'ibs_individual_matrix_260522.csv'))

print("\nAll 21 genetic analyses completed successfully and saved to `/output/rawvisual/Tib/`!")
