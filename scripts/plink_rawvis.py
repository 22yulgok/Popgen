"""
PLINK Genotype Exploratory Visualization
=========================================
Input : /raw/Tibetan.fstat.1240K.260508.{bed,bim,fam}
Output: /output/rawvis/  (12 PNG figures + optional CSV summaries)

BED parsing: manual numpy-based SNP-major decoding (no PLINK binary required).
Genotype encoding (SNP-major):
  2-bit code | meaning
  ---------- | --------
  00 (0)     | hom-ref  (A1/A1)
  01 (1)     | MISSING
  10 (2)     | heterozygous (A1/A2)
  11 (3)     | hom-alt  (A2/A2)
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from collections import Counter

warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────────────
PREFIX   = '/home/duri_bae/popgen_lab/raw/Tibetan.fstat.1240K.260508'
FAM_PATH = PREFIX + '.fam'
BIM_PATH = PREFIX + '.bim'
BED_PATH = PREFIX + '.bed'
OUT_DIR  = '/home/duri_bae/popgen_lab/output/rawvis'
OUT_DIR_ALT = '/home/duri_bae/popgen_lab/output/plinkvisual'
DS       = '260522'   # date stamp

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(OUT_DIR_ALT, exist_ok=True)

def save(name):
    p = os.path.join(OUT_DIR, f'{name}_{DS}.png')
    plt.savefig(p, dpi=150, bbox_inches='tight')
    
    # Save to the alternate directory too
    import shutil
    p_alt = os.path.join(OUT_DIR_ALT, f'{name}_{DS}.png')
    shutil.copy(p, p_alt)
    
    plt.close()
    print(f'  → {p}')
    print(f'  → {p_alt}')
    return p

plt.rcParams.update({
    'figure.dpi': 150,
    'font.family': 'DejaVu Sans',
    'font.size': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.titlesize': 12,
    'axes.titleweight': 'bold',
    'axes.labelsize': 10,
    'legend.fontsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
})

# ══════════════════════════════════════════════════════════════════════════════
# 1 — Load FAM and BIM (lightweight, no BED needed)
# ══════════════════════════════════════════════════════════════════════════════
print('Loading FAM and BIM …')
fam = pd.read_csv(FAM_PATH, sep=r'\s+', header=None,
                  names=['FID','IID','PID','MID','SEX','PHENO'])
bim = pd.read_csv(BIM_PATH, sep=r'\t', header=None,
                  names=['CHR','SNP','CM','POS','A1','A2'])

N_SAM = len(fam)
N_SNP = len(bim)
BYTES_PER_SNP = -(-N_SAM // 4)   # ceil(N_SAM / 4)

print(f'  Samples : {N_SAM}')
print(f'  SNPs    : {N_SNP:,}')
print(f'  FIDs    : {fam.FID.nunique()}')
print(f'  Chroms  : {sorted(bim.CHR.unique())}')

# ══════════════════════════════════════════════════════════════════════════════
# 2 — Parse BED: compute per-sample and per-SNP missingness + alt allele count
#     Process in chunks of CHUNK_SIZE SNPs to avoid RAM overflow
# ══════════════════════════════════════════════════════════════════════════════
CHUNK_SIZE = 5000

print(f'\nParsing BED in chunks of {CHUNK_SIZE} SNPs …')

# Precompute lookup table: byte → 4 genotype codes (0–3)
lut = np.zeros((256, 4), dtype=np.uint8)
for byte in range(256):
    for bit_pair in range(4):
        lut[byte, bit_pair] = (byte >> (2 * bit_pair)) & 3

# Accumulators
sample_missing   = np.zeros(N_SAM, dtype=np.int32)   # missing calls per sample
sample_allele_ct = np.zeros(N_SAM, dtype=np.int64)   # total called alleles (non-missing, diploid)
sample_alt_ct    = np.zeros(N_SAM, dtype=np.int64)   # alt allele count per sample (for sanity)

snp_missing  = np.zeros(N_SNP, dtype=np.int32)       # missing calls per SNP
snp_alt_hom  = np.zeros(N_SNP, dtype=np.int32)       # count of hom-alt (11)
snp_het      = np.zeros(N_SNP, dtype=np.int32)       # count of het (10)
snp_ref_hom  = np.zeros(N_SNP, dtype=np.int32)       # count of hom-ref (00)

with open(BED_PATH, 'rb') as f:
    magic = f.read(3)
    assert magic == b'\x6c\x1b\x01', f'Unexpected BED magic: {magic.hex()}'

    snp_idx = 0
    while snp_idx < N_SNP:
        chunk_end = min(snp_idx + CHUNK_SIZE, N_SNP)
        chunk_n   = chunk_end - snp_idx

        # Read raw bytes for this chunk
        raw = np.frombuffer(f.read(chunk_n * BYTES_PER_SNP), dtype=np.uint8)
        raw = raw.reshape(chunk_n, BYTES_PER_SNP)

        # Decode: shape (chunk_n, BYTES_PER_SNP, 4) → (chunk_n, BYTES_PER_SNP*4)
        decoded = lut[raw]                               # (chunk_n, bytes, 4)
        decoded = decoded.reshape(chunk_n, -1)[:, :N_SAM]  # (chunk_n, N_SAM)

        missing_mask = (decoded == 1)   # MISSING
        ref_hom_mask = (decoded == 0)   # 00
        het_mask     = (decoded == 2)   # 10
        alt_hom_mask = (decoded == 3)   # 11

        # Per-sample missingness accumulation
        sample_missing += missing_mask.sum(axis=0)

        # Per-SNP stats
        snp_missing[snp_idx:chunk_end]  = missing_mask.sum(axis=1)
        snp_ref_hom[snp_idx:chunk_end]  = ref_hom_mask.sum(axis=1)
        snp_het[snp_idx:chunk_end]      = het_mask.sum(axis=1)
        snp_alt_hom[snp_idx:chunk_end]  = alt_hom_mask.sum(axis=1)

        snp_idx = chunk_end
        if snp_idx % 100_000 == 0:
            print(f'    parsed {snp_idx:,} / {N_SNP:,} SNPs …')

print(f'  BED parsing complete.')

# Derived statistics
sample_callrate = 1.0 - sample_missing / N_SNP
snp_callrate    = 1.0 - snp_missing / N_SAM

# MAF from snp_ref_hom, snp_het, snp_alt_hom
# A2 allele count = 2*snp_alt_hom + snp_het
# Total allele count = 2*(N_SAM - snp_missing)
called_n  = N_SAM - snp_missing.astype(np.float64)
total_al  = 2 * called_n
a2_ct     = 2 * snp_alt_hom.astype(np.float64) + snp_het.astype(np.float64)
a2_freq   = np.where(total_al > 0, a2_ct / total_al, np.nan)
maf       = np.where(a2_freq > 0.5, 1.0 - a2_freq, a2_freq)

print(f'  Sample call rate: mean={sample_callrate.mean():.4f}  min={sample_callrate.min():.4f}')
print(f'  SNP call rate:    mean={snp_callrate.mean():.4f}  min={snp_callrate.min():.4f}')
print(f'  MAF:              mean={np.nanmean(maf):.4f}  monomorphic={np.sum(maf == 0):,}')

# ── Attach SNP stats back to bim ─────────────────────────────────────────────
bim['missing']   = snp_missing
bim['callrate']  = snp_callrate
bim['maf']       = maf
bim['ref_hom']   = snp_ref_hom
bim['het']       = snp_het
bim['alt_hom']   = snp_alt_hom

fam['missing']   = sample_missing
fam['callrate']  = sample_callrate

# ── Save summary CSVs ─────────────────────────────────────────────────────────
bim_out = bim[['CHR','SNP','POS','callrate','maf']].copy()
bim_out.to_csv(os.path.join(OUT_DIR, f'variant_summary_{DS}.csv'), index=False)
bim_out.to_csv(os.path.join(OUT_DIR_ALT, f'variant_summary_{DS}.csv'), index=False)

fam_out = fam[['FID','IID','SEX','PHENO','callrate']].copy()
fam_out.to_csv(os.path.join(OUT_DIR, f'sample_summary_{DS}.csv'), index=False)
fam_out.to_csv(os.path.join(OUT_DIR_ALT, f'sample_summary_{DS}.csv'), index=False)

print(f'\nSummary CSVs saved to {OUT_DIR} and {OUT_DIR_ALT}')

# ══════════════════════════════════════════════════════════════════════════════
# PALETTE
# ══════════════════════════════════════════════════════════════════════════════
CHR_COLORS = plt.cm.tab20.colors   # 20 colours, reused

# ══════════════════════════════════════════════════════════════════════════════
# Figure 1 — SNP count by chromosome
# ══════════════════════════════════════════════════════════════════════════════
print('\n[1/12] SNP count by chromosome …')

chr_counts = bim['CHR'].value_counts().sort_index()
chr_labels = {23: 'X', 24: 'Y'}
labels = [chr_labels.get(c, str(c)) for c in chr_counts.index]

fig, ax = plt.subplots(figsize=(12, 5))
bars = ax.bar(labels, chr_counts.values,
              color=[CHR_COLORS[i % 20] for i in range(len(chr_counts))],
              edgecolor='none')
ax.set_xlabel('Chromosome')
ax.set_ylabel('SNP count')
ax.set_title(f'SNP Count by Chromosome\n(total {N_SNP:,} SNPs, {N_SAM} samples)')
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{int(x):,}'))
for bar, val in zip(bars, chr_counts.values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 200,
            f'{val:,}', ha='center', va='bottom', fontsize=6.5)
plt.tight_layout()
save('snp_count_by_chr')

# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 — SNP density (inter-SNP spacing) by chromosome
# ══════════════════════════════════════════════════════════════════════════════
print('[2/12] SNP density by chromosome …')

# Subsample for speed: use every 10th SNP per chr to compute spacing
spacing_data = {}
for chrom, grp in bim.groupby('CHR'):
    pos_sorted = np.sort(grp['POS'].values)
    if len(pos_sorted) > 1:
        diffs = np.diff(pos_sorted[::10])   # every 10th position
        spacing_data[chr_labels.get(chrom, str(chrom))] = diffs

chr_order = [chr_labels.get(c, str(c)) for c in sorted(bim['CHR'].unique())]
sp_filtered = [spacing_data[c] for c in chr_order if c in spacing_data]

fig, ax = plt.subplots(figsize=(13, 5))
bp = ax.boxplot(sp_filtered, labels=chr_order,
                patch_artist=True, notch=False, showfliers=False,
                medianprops=dict(color='black', lw=1.5))
for patch, color in zip(bp['boxes'], [CHR_COLORS[i % 20] for i in range(len(chr_order))]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_xlabel('Chromosome')
ax.set_ylabel('Inter-SNP spacing (bp, log scale)')
ax.set_yscale('log')
ax.set_title('SNP Density: Inter-SNP Spacing by Chromosome\n(sampled every 10th SNP per chromosome, outliers hidden)')
plt.tight_layout()
save('snp_density_by_chr')

# ══════════════════════════════════════════════════════════════════════════════
# Figure 3 — Physical position distribution (Manhattan-style)
# ══════════════════════════════════════════════════════════════════════════════
print('[3/12] Physical position distribution …')

# Subsample 2% for readability
rng = np.random.default_rng(42)
sub_idx = rng.choice(N_SNP, size=N_SNP // 50, replace=False)
sub_bim = bim.iloc[sub_idx].copy()

fig, ax = plt.subplots(figsize=(14, 4))
x_offset = 0
xtick_pos, xtick_lab = [], []
for i, (chrom, grp) in enumerate(sub_bim.groupby('CHR')):
    x = grp['POS'].values / 1e6 + x_offset
    ax.scatter(x, [i % 2 * 0.15 + 0.4] * len(x),
               s=0.3, color=CHR_COLORS[i % 20], alpha=0.4, rasterized=True)
    xtick_pos.append(x_offset + grp['POS'].values.max() / 1e6 / 2)
    xtick_lab.append(chr_labels.get(chrom, str(chrom)))
    x_offset += grp['POS'].values.max() / 1e6 + 5

ax.set_xticks(xtick_pos)
ax.set_xticklabels(xtick_lab, fontsize=7)
ax.set_yticks([])
ax.set_xlabel('Genomic position (cumulative, Mb)')
ax.set_title(f'Variant Position Distribution (2% subsample, n={len(sub_bim):,})\nColoured by chromosome')
ax.set_xlim(-5, x_offset)
plt.tight_layout()
save('position_distribution')

# ══════════════════════════════════════════════════════════════════════════════
# Figure 4 — Allele type distribution
# ══════════════════════════════════════════════════════════════════════════════
print('[4/12] Allele type distribution …')

bim['allele_pair'] = bim.apply(
    lambda r: '/'.join(sorted([str(r['A1']).upper(), str(r['A2']).upper()])), axis=1)
pair_counts = bim['allele_pair'].value_counts()

fig, ax = plt.subplots(figsize=(8, 5))
colors_pair = plt.cm.Set2.colors[:len(pair_counts)]
bars = ax.bar(pair_counts.index, pair_counts.values,
              color=colors_pair, edgecolor='none')
for bar, val in zip(bars, pair_counts.values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 500,
            f'{val:,}\n({val/N_SNP*100:.1f}%)', ha='center', va='bottom', fontsize=8)
ax.set_xlabel('Allele pair (sorted)')
ax.set_ylabel('SNP count')
ax.set_title(f'Allele Type Distribution\n({N_SNP:,} total SNPs)')
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{int(x):,}'))
plt.tight_layout()
save('allele_type_distribution')

# ══════════════════════════════════════════════════════════════════════════════
# Figure 5 — Transition / Transversion ratio by chromosome
# ══════════════════════════════════════════════════════════════════════════════
print('[5/12] Ti/Tv ratio by chromosome …')

TRANSITIONS = {frozenset({'A','G'}), frozenset({'C','T'})}
def is_ti(a1, a2):
    return frozenset({str(a1).upper(), str(a2).upper()}) in TRANSITIONS

bim['is_ti'] = bim.apply(lambda r: is_ti(r['A1'], r['A2']), axis=1)

titv_chr = bim.groupby('CHR')['is_ti'].agg(
    ti='sum', tv=lambda x: (~x).sum())
titv_chr['ratio'] = titv_chr['ti'] / titv_chr['tv'].replace(0, np.nan)
titv_chr.index = [chr_labels.get(c, str(c)) for c in titv_chr.index]

fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
x = range(len(titv_chr))

axes[0].bar(x, titv_chr['ti'], label='Transitions', color='#4393C3', alpha=0.8)
axes[0].bar(x, titv_chr['tv'], bottom=titv_chr['ti'], label='Transversions',
            color='#D73027', alpha=0.8)
axes[0].set_ylabel('SNP count')
axes[0].set_title('Ti/Tv Counts and Ratio by Chromosome')
axes[0].legend()
axes[0].yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{int(x):,}'))

axes[1].bar(x, titv_chr['ratio'], color='#4DAC26', alpha=0.8)
axes[1].axhline(titv_chr['ratio'].mean(), color='black', lw=1.5, ls='--',
                label=f'Mean Ti/Tv = {titv_chr["ratio"].mean():.3f}')
axes[1].set_ylabel('Ti/Tv ratio')
axes[1].set_xticks(list(x))
axes[1].set_xticklabels(titv_chr.index, fontsize=8)
axes[1].set_xlabel('Chromosome')
axes[1].legend()

plt.tight_layout()
save('ti_tv_ratio')

# ══════════════════════════════════════════════════════════════════════════════
# Figure 6 — Sex code distribution
# ══════════════════════════════════════════════════════════════════════════════
print('[6/12] Sex code distribution …')

sex_map = {0: 'Unknown', 1: 'Male', 2: 'Female'}
sex_counts = fam['SEX'].map(sex_map).value_counts()

fig, ax = plt.subplots(figsize=(6, 6))
colors_sex = {'Male': '#2166AC', 'Female': '#D73027', 'Unknown': '#888888'}
wedge_cols = [colors_sex.get(k, '#aaa') for k in sex_counts.index]
wedges, texts, autotexts = ax.pie(
    sex_counts.values, labels=sex_counts.index,
    autopct=lambda p: f'{p:.1f}%\n(n={int(round(p*N_SAM/100))})',
    colors=wedge_cols, startangle=90,
    textprops={'fontsize': 11})
ax.set_title(f'Sex Code Distribution\n(n={N_SAM} samples)')
plt.tight_layout()
save('sex_code_distribution')

# ══════════════════════════════════════════════════════════════════════════════
# Figure 7 — Sample count by FID (population-like group)
# ══════════════════════════════════════════════════════════════════════════════
print('[7/12] Sample count by population (FID) …')

pop_counts = fam['FID'].value_counts().sort_values()

fig, ax = plt.subplots(figsize=(10, max(12, len(pop_counts) * 0.22)))
bar_colors = [CHR_COLORS[i % 20] for i in range(len(pop_counts))]
ax.barh(pop_counts.index, pop_counts.values, color=bar_colors, edgecolor='none', height=0.75)
ax.set_xlabel('Sample count')
ax.set_title(f'Sample Count by FID Group\n(FID treated as population-like label, n={N_SAM} samples, {fam.FID.nunique()} groups)',
             fontweight='bold')
ax.axvline(pop_counts.mean(), color='black', lw=1.2, ls='--',
           label=f'Mean = {pop_counts.mean():.1f}')
ax.legend(fontsize=8)
for val, name in zip(pop_counts.values, pop_counts.index):
    ax.text(val + 0.2, name, str(val), va='center', fontsize=6)
plt.tight_layout()
save('sample_count_by_pop')

# ══════════════════════════════════════════════════════════════════════════════
# Figure 8 — Per-sample missingness distribution
# ══════════════════════════════════════════════════════════════════════════════
print('[8/12] Per-sample missingness distribution …')

miss_rate = 1 - sample_callrate

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
ax.hist(miss_rate, bins=60, color='#4393C3', edgecolor='white', lw=0.3)
ax.axvline(miss_rate.mean(), color='#D73027', lw=1.5, ls='--',
           label=f'Mean = {miss_rate.mean():.4f}')
ax.axvline(0.10, color='#F4A582', lw=1.2, ls=':', label='10% threshold')
ax.set_xlabel('Missing rate')
ax.set_ylabel('Sample count')
ax.set_title('Per-sample Missingness Distribution')
ax.legend()

ax2 = axes[1]
ax2.hist(sample_callrate, bins=60, color='#4393C3', edgecolor='white', lw=0.3)
ax2.axvline(sample_callrate.mean(), color='#D73027', lw=1.5, ls='--',
            label=f'Mean call rate = {sample_callrate.mean():.4f}')
ax2.set_xlabel('Call rate')
ax2.set_ylabel('Sample count')
ax2.set_title('Per-sample Call Rate Distribution')
ax2.legend()

n_fail_miss = int((miss_rate > 0.10).sum())
fig.suptitle(f'Per-sample Missingness (n={N_SAM} samples)\n'
             f'{n_fail_miss} samples >10% missing',
             fontweight='bold')
plt.tight_layout()
save('sample_missingness_dist')

# ══════════════════════════════════════════════════════════════════════════════
# Figure 9 — Per-variant missingness distribution
# ══════════════════════════════════════════════════════════════════════════════
print('[9/12] Per-variant missingness distribution …')

var_miss = 1 - snp_callrate

fig, ax = plt.subplots(figsize=(10, 5))
ax.hist(var_miss, bins=100, color='#74ADD1', edgecolor='white', lw=0.2, log=True)
ax.axvline(var_miss.mean(), color='#D73027', lw=1.5, ls='--',
           label=f'Mean = {var_miss.mean():.4f}')
ax.axvline(0.05, color='#F4A582', lw=1.2, ls=':', label='5% threshold')
ax.axvline(0.10, color='#FDAE61', lw=1.2, ls=':', label='10% threshold')
n_fail5 = int((var_miss > 0.05).sum())
n_fail10 = int((var_miss > 0.10).sum())
ax.set_xlabel('Missing rate per variant')
ax.set_ylabel('Variant count (log scale)')
ax.set_title(f'Per-variant Missingness Distribution\n'
             f'{N_SNP:,} total variants | >5%: {n_fail5:,} | >10%: {n_fail10:,}')
ax.legend()
plt.tight_layout()
save('variant_missingness_dist')

# ══════════════════════════════════════════════════════════════════════════════
# Figure 10 — Variant missingness by chromosome (box plot)
# ══════════════════════════════════════════════════════════════════════════════
print('[10/12] Variant missingness by chromosome …')

miss_by_chr = [bim[bim['CHR'] == c]['callrate'].values
               for c in sorted(bim['CHR'].unique())]
chr_lab_ordered = [chr_labels.get(c, str(c)) for c in sorted(bim['CHR'].unique())]

fig, ax = plt.subplots(figsize=(13, 5))
bp = ax.boxplot(miss_by_chr, labels=chr_lab_ordered,
                patch_artist=True, notch=False, showfliers=False,
                medianprops=dict(color='black', lw=1.5))
for patch, color in zip(bp['boxes'], [CHR_COLORS[i % 20] for i in range(len(chr_lab_ordered))]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.axhline(snp_callrate.mean(), color='#D73027', lw=1.2, ls='--',
           label=f'Overall mean call rate = {snp_callrate.mean():.4f}')
ax.set_xlabel('Chromosome')
ax.set_ylabel('Variant call rate')
ax.set_title('Variant Call Rate Distribution by Chromosome\n(outliers hidden)')
ax.legend()
plt.tight_layout()
save('missingness_by_chr')

# ══════════════════════════════════════════════════════════════════════════════
# Figure 11 — MAF distribution
# ══════════════════════════════════════════════════════════════════════════════
print('[11/12] MAF distribution …')

maf_valid = maf[~np.isnan(maf)]
n_mono     = int((maf_valid == 0).sum())
n_rare     = int(((maf_valid > 0) & (maf_valid < 0.05)).sum())
n_common   = int((maf_valid >= 0.05).sum())

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle(f'Minor Allele Frequency Distribution\n'
             f'{N_SNP:,} total SNPs | monomorphic: {n_mono:,} | MAF<5%: {n_rare:,} | MAF≥5%: {n_common:,}',
             fontweight='bold')

ax = axes[0]
ax.hist(maf_valid, bins=100, color='#4393C3', edgecolor='white', lw=0.2)
ax.axvline(0.05, color='#F4A582', lw=1.2, ls='--', label='MAF = 0.05')
ax.set_xlabel('MAF')
ax.set_ylabel('SNP count')
ax.set_title('Full MAF distribution')
ax.legend()

ax2 = axes[1]
maf_nonzero = maf_valid[maf_valid > 0]
ax2.hist(maf_nonzero, bins=100, color='#4393C3', edgecolor='white', lw=0.2)
ax2.axvline(0.05, color='#F4A582', lw=1.2, ls='--', label='MAF = 0.05')
ax2.axvline(maf_nonzero.mean(), color='#D73027', lw=1.2, ls='--',
            label=f'Mean MAF = {maf_nonzero.mean():.4f}')
ax2.set_xlabel('MAF')
ax2.set_ylabel('SNP count')
ax2.set_title('MAF > 0 (excluding monomorphic)')
ax2.legend()

plt.tight_layout()
save('maf_distribution')

# ══════════════════════════════════════════════════════════════════════════════
# Figure 12 — Top missing samples (ranked bar)
# ══════════════════════════════════════════════════════════════════════════════
print('[12/12] Top missing samples …')

fam_sorted = fam.sort_values('callrate').head(40).copy()
fam_sorted['label'] = fam_sorted['FID'] + '/' + fam_sorted['IID']
miss_pct = (1 - fam_sorted['callrate']) * 100

fig, ax = plt.subplots(figsize=(10, 10))
bar_colors = ['#D73027' if v > 10 else '#4393C3' for v in miss_pct.values]
ax.barh(fam_sorted['label'], miss_pct.values, color=bar_colors, edgecolor='none', height=0.75)
ax.axvline(10, color='#F4A582', lw=1.2, ls='--', label='>10% missing threshold')
ax.set_xlabel('Missing rate (%)')
ax.set_title('Top 40 Samples by Missingness\n(red = >10% missing)')
ax.legend()
for val, lab in zip(miss_pct.values, fam_sorted['label']):
    ax.text(val + 0.05, lab, f'{val:.1f}%', va='center', fontsize=7)
plt.tight_layout()
save('top_missing_samples')

# ══════════════════════════════════════════════════════════════════════════════
# FINAL REPORT
# ══════════════════════════════════════════════════════════════════════════════
generated = sorted(f for f in os.listdir(OUT_DIR) if f.endswith('.png'))
csvs      = sorted(f for f in os.listdir(OUT_DIR) if f.endswith('.csv'))

print('\n' + '='*65)
print('FINAL REPORT')
print('='*65)
print(f'\nInput files:')
print(f'  FAM: {FAM_PATH}  — OK ({N_SAM} samples)')
print(f'  BIM: {BIM_PATH}  — OK ({N_SNP:,} SNPs)')
print(f'  BED: {BED_PATH}  — OK (magic 6c1b01, SNP-major)')
print(f'\nDataset summary:')
print(f'  Samples           : {N_SAM}')
print(f'  SNPs              : {N_SNP:,}')
print(f'  Chromosomes       : 1–22, 23 (X), 24 (Y)')
print(f'  Unique FIDs       : {fam.FID.nunique()}  (treated as pop-like labels)')
print(f'  Sex (M/F/unknown) : {(fam.SEX==1).sum()} / {(fam.SEX==2).sum()} / {(fam.SEX==0).sum()}')
print(f'  Phenotype         : all coded 2 (uninformative)')
print(f'  Mean sample call rate : {sample_callrate.mean():.5f}')
print(f'  Mean SNP call rate    : {snp_callrate.mean():.5f}')
print(f'  Mean MAF (non-mono)   : {np.nanmean(maf[maf > 0]):.5f}')
print(f'  Monomorphic SNPs      : {n_mono:,}')
print(f'\nGenerated {len(generated)} figures:')
for f in generated:
    print(f'  {f}')
print(f'\nGenerated {len(csvs)} CSV summaries:')
for f in csvs:
    print(f'  {f}')
print(f'\nAssumptions:')
print(f'  • FID used as population-like label (not necessarily true population ID)')
print(f'  • Phenotype code 2 treated as uninformative (not disease status)')
print(f'  • Chr 23 = X chromosome, chr 24 = Y chromosome (standard PLINK coding)')
print(f'  • Position distribution subsample: 2% random (for readability)')
print(f'  • SNP density: every 10th position sampled per chromosome')
print(f'\nSkipped plots:')
print(f'  • PCA: computationally heavy for 1M+ SNPs without PLINK; skipped.')
print(f'    Alternative: use PLINK --pca with LD pruning first.')
print(f'  • LD pruning summary: PLINK not available; skipped.')
print('\nDone.')
