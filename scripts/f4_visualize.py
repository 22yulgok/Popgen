"""
F4 Visualization Script — Tibetan Ancient DNA Analysis
=======================================================
Resolves dates for ALL Pop2 groups (including unlabelled ones like SDLG, GBSL,
aMMD_*, MBC*) by matching IIDs from the .ind file into the AADR "Genetic ID"
column, then averaging per PID group.

nSNPs cutoff: >= 100,000 (applied before any plotting)
Output: /output/f4visual/
"""

import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
import seaborn as sns
from scipy.stats import gaussian_kde

warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────────────
F4_PATH   = '/home/duri_bae/popgen_lab/output/Tibetan.f4.1240K.260506.csv'
IND_PATH  = '/home/duri_bae/popgen_lab/raw/Tibetan.fstat.1240K.260508.ind'
SI_PATH   = '/home/duri_bae/popgen_lab/raw/SampleInfo_260522.csv'
AADR_PATH = '/home/duri_bae/popgen_lab/raw/v66.1240K.aadr.PUB.csv'
OUT_DIR   = '/home/duri_bae/popgen_lab/output/f4visual'
DATE_STAMP = '20260522'
NSNP_CUTOFF = 100_000

os.makedirs(OUT_DIR, exist_ok=True)

# ── Style ──────────────────────────────────────────────────────────────────────
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
# Step 1 — Load and filter the f4 result table
# ══════════════════════════════════════════════════════════════════════════════
print("Loading f4 data …")
df = pd.read_csv(F4_PATH)
df.columns = df.columns.str.strip()
df = df.rename(columns={'f4.est': 'f4est'})
n_raw = len(df)

df = df[df['nSNPs'] >= NSNP_CUTOFF].copy()
df_noNA = df.dropna(subset=['z']).copy()
print(f"  Raw rows:       {n_raw:,}")
print(f"  After nSNPs≥{NSNP_CUTOFF//1000}K: {len(df):,}  (removed {n_raw - len(df):,})")
print(f"  After z NA drop: {len(df_noNA):,}")

# ══════════════════════════════════════════════════════════════════════════════
# Step 2 — Build a PID → mean_date_BP map using .ind + AADR
# ══════════════════════════════════════════════════════════════════════════════
print("\nBuilding PID → date map from AADR …")
ind = pd.read_csv(IND_PATH, sep=r'\s+', header=None, names=['IID', 'SEX', 'PID'])

aadr = pd.read_csv(AADR_PATH, dtype=str, low_memory=False)
genetic_id_col = aadr.columns[0]   # Genetic ID (with suffixes)
group_id_col   = aadr.columns[14]  # Group ID  ≈ PID
date_col       = aadr.columns[10]  # Date mean in BP
aadr['Date_Numeric'] = pd.to_numeric(aadr[date_col], errors='coerce')

# Merge ind → AADR on Genetic ID
ind_aadr = ind.merge(
    aadr[[genetic_id_col, group_id_col, 'Date_Numeric']],
    left_on='IID', right_on=genetic_id_col, how='left'
)

# Average date per PID (using only rows where date was found)
pid_date_map = (
    ind_aadr.dropna(subset=['Date_Numeric'])
    .groupby('PID')['Date_Numeric']
    .mean()
    .to_dict()
)
print(f"  PIDs with date resolved: {len(pid_date_map)}")

# Also load SampleInfo for period labels (PID may not be unique — take first)
si = pd.read_csv(SI_PATH)
si_pid_period = si.groupby('PID')['Period'].first().to_dict()
si_pid_cal    = si.groupby('PID')['Cal'].first().to_dict()

# ══════════════════════════════════════════════════════════════════════════════
# Step 3 — Annotate all Pop2 groups with date and geographic region
# ══════════════════════════════════════════════════════════════════════════════

# Known regional classification (by population name prefix / research context)
TP_PREFIXES = [
    'Zongri', 'Yushu', 'Shannan', 'Chamdo', 'Shigatse', 'Lhasa',
    'Nagqu', 'Ngari', 'Nyingchi', 'SDLG', 'GBSL', 'QLSZ', 'PLDW', 'Guge', 'Laga',
]
HIMALAYAN_POPS = ['Chokhopani', 'Rhirhi', 'Kyang', 'Mebrak', 'Samdzong', 'Lubrak', 'Suila']
MBC_PREFIXES   = ['Mbc', 'Mabuco', 'MabucoE', 'Xingyi', 'DHL_', 'Haneyi', 'BGD_EN']

def classify_pop(pop):
    for pfx in TP_PREFIXES:
        if pop.startswith(pfx):
            return 'Tibetan_Plateau'
    if pop in HIMALAYAN_POPS:
        return 'Himalayan'
    for pfx in MBC_PREFIXES:
        if pop.startswith(pfx):
            return 'Sichuan_Mabuco'
    return 'Other'

all_pop2 = df['Pop2'].unique()
pop2_meta = {}
for pop in all_pop2:
    region = classify_pop(pop)
    # Use .ind-based date lookup
    date_bp = pid_date_map.get(pop, None)
    # Fallback: extract trailing numeric date from label, e.g. Zongri5.1k → 5100
    if date_bp is None:
        m = re.search(r'(\d+\.?\d*)k', pop, re.IGNORECASE)
        if m:
            date_bp = round(float(m.group(1)) * 1000)
    pop2_meta[pop] = {'region': region, 'date_bp': date_bp}

# Build tidy DataFrame of pop2 metadata
meta_df = pd.DataFrame.from_dict(pop2_meta, orient='index').reset_index()
meta_df.columns = ['Pop2', 'region', 'date_bp']
matched_count   = meta_df['date_bp'].notna().sum()
unmatched_count = meta_df['date_bp'].isna().sum()
print(f"  Pop2 groups with date resolved:   {matched_count}")
print(f"  Pop2 groups with date unknown:    {unmatched_count}  (metadata_not_found)")

df = df.merge(meta_df, on='Pop2', how='left')
df_noNA = df.dropna(subset=['z']).copy()

# Shorthand selectors
tp_pops   = meta_df[meta_df['region'] == 'Tibetan_Plateau']['Pop2'].tolist()
him_pops  = meta_df[meta_df['region'] == 'Himalayan']['Pop2'].tolist()
mbc_pops  = meta_df[meta_df['region'] == 'Sichuan_Mabuco']['Pop2'].tolist()

# ══════════════════════════════════════════════════════════════════════════════
# REGION COLOR MAP
# ══════════════════════════════════════════════════════════════════════════════
REGION_COLORS = {
    'Tibetan_Plateau': '#2166AC',
    'Himalayan'      : '#74ADD1',
    'Sichuan_Mabuco' : '#D73027',
    'Other'          : '#888888',
}

# ══════════════════════════════════════════════════════════════════════════════
# Figure 1: Z-score distribution (all filtered results)
# ══════════════════════════════════════════════════════════════════════════════
print("\n[1/7] Z-score distribution …")

z = df_noNA['z']
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle(
    f'f4 Z-score Distribution (nSNPs ≥ {NSNP_CUTOFF//1000}K, n={len(z):,})\n'
    'f4(Mbuti.DG, Pop2; Pop3, Pop4)',
    fontsize=12, fontweight='bold'
)

# Full distribution
ax = axes[0]
ax.hist(z, bins=200, color='#4393C3', edgecolor='none', alpha=0.7, density=True)
kde_x = np.linspace(z.min(), z.max(), 1000)
kde_y = gaussian_kde(z.sample(min(30000, len(z)), random_state=42))(kde_x)
ax.plot(kde_x, kde_y, color='#D73027', lw=1.5, label='KDE')
for thr, col in [(-3, '#F4A582'), (3, '#F4A582'), (0, 'black')]:
    ax.axvline(thr, color=col, lw=1.2 if thr == 0 else 1.5,
               ls='-' if thr == 0 else '--',
               alpha=0.5 if thr == 0 else 1.0,
               label='|Z| = 3 threshold' if thr == 3 else None)
sig_pct = (np.abs(z) >= 3).mean() * 100
ax.set_xlabel('Z-score'); ax.set_ylabel('Density')
ax.set_title(f'Full range\n{sig_pct:.1f}% have |Z| \u2265 3')
ax.legend(loc='upper right')

# Zoom ±15
ax2 = axes[1]
z_c = z.clip(-15, 15)
ax2.hist(z_c, bins=150, color='#4393C3', edgecolor='none', alpha=0.7)
ax2.axvline(-3, color='#F4A582', lw=1.5, ls='--', label='|Z| = 3')
ax2.axvline( 3, color='#F4A582', lw=1.5, ls='--')
ax2.axvline( 0, color='black', lw=1.0, alpha=0.5)
ax2.set_xlabel('Z-score (clipped ±15)'); ax2.set_ylabel('Count')
ax2.set_title('Zoom: Z in [−15, 15]'); ax2.legend(loc='upper right')

plt.tight_layout()
p = os.path.join(OUT_DIR, f'f4_zscore_distribution_{DATE_STAMP}.png')
plt.savefig(p, bbox_inches='tight'); plt.close(); print(f"  → {p}")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 2: Top 40 most extreme Z-scores among Tibetan Plateau Pop2
# ══════════════════════════════════════════════════════════════════════════════
print("[2/7] Ranked significant results …")

tp_sig = df_noNA[df_noNA['Pop2'].isin(tp_pops) & (np.abs(df_noNA['z']) >= 3)].copy()
top40 = pd.concat([
    tp_sig.nlargest(20, 'z'),
    tp_sig.nsmallest(20, 'z'),
]).drop_duplicates().sort_values('z')

fig, ax = plt.subplots(figsize=(11, 12))
colors = ['#313695' if v < 0 else '#D73027' for v in top40['z']]
ax.barh(range(len(top40)), top40['z'], color=colors, edgecolor='none', height=0.72)
ax.axvline(0, color='black', lw=0.8)
ax.axvline(-3, color='#F4A582', lw=1.2, ls='--', alpha=0.8, label='|Z| = 3')
ax.axvline( 3, color='#F4A582', lw=1.2, ls='--', alpha=0.8)
labels = [f"{r['Pop2']}  |  {r['Pop3']} vs {r['Pop4']}" for _, r in top40.iterrows()]
ax.set_yticks(range(len(top40))); ax.set_yticklabels(labels, fontsize=7.5)
ax.set_xlabel('Z-score')
ax.set_title(
    f'Top 40 Extreme f4 Results — Tibetan Plateau Pop2 (|Z| ≥ 3, nSNPs ≥ {NSNP_CUTOFF//1000}K)\n'
    'f4(Mbuti.DG, Pop2; Pop3, Pop4)',
    fontweight='bold'
)
pos_p = mpatches.Patch(color='#D73027', label='Positive Z (closer to Pop3)')
neg_p = mpatches.Patch(color='#313695', label='Negative Z (closer to Pop4)')
thr_l = plt.Line2D([0],[0], color='#F4A582', ls='--', lw=1.2, label='|Z| = 3')
ax.legend(handles=[pos_p, neg_p, thr_l], loc='lower right', fontsize=8)
plt.tight_layout()
p = os.path.join(OUT_DIR, f'f4_significant_results_{DATE_STAMP}.png')
plt.savefig(p, bbox_inches='tight'); plt.close(); print(f"  → {p}")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 3: Forest plot — Tibetan Plateau sites vs (Xingyi_EN, YR_MN)
# ══════════════════════════════════════════════════════════════════════════════
print("[3/7] Forest plot: TP sites vs Xingyi_EN / YR_MN …")

P3, P4 = 'Xingyi_EN', 'YR_MN'
sub3 = df[(df['Pop3'] == P3) & (df['Pop4'] == P4) & df['Pop2'].isin(tp_pops)].copy()
# Sort: date descending (older on top), unknown date at bottom
sub3['date_sort'] = sub3['date_bp'].fillna(-1)
sub3 = sub3.sort_values('date_sort', ascending=False)

# Sub-region colour by name prefix
def tp_subregion(pop):
    if pop.startswith('Zongri'):    return '#1A6FBF'   # NE Tibetan
    if pop.startswith('Yushu'):     return '#6A3D9A'
    if pop.startswith('Shannan'):   return '#33A02C'
    if pop.startswith('Chamdo'):    return '#E31A1C'
    if pop.startswith('Shigatse'):  return '#FF7F00'
    if pop.startswith('Lhasa'):     return '#F4A500'
    if pop.startswith('Nagqu'):     return '#B15928'
    if pop.startswith('Ngari'):     return '#1F78B4'
    if pop.startswith('Nyingchi'):  return '#A6CEE3'
    return '#888888'   # modern/other TP

fig, ax = plt.subplots(figsize=(9, max(7, len(sub3) * 0.33)))
for i, (_, row) in enumerate(sub3.iterrows()):
    c = tp_subregion(row['Pop2'])
    ax.errorbar(row['f4est'], i, xerr=1.96*row['se'],
                fmt='o', color=c, ecolor=c, ms=5, elinewidth=1.3, capsize=2.5)

ax.axvline(0, color='black', lw=0.8)
ax.axvline(-3*sub3['se'].median(), color='#F4A582', lw=1, ls='--', alpha=0.5)
ax.axvline( 3*sub3['se'].median(), color='#F4A582', lw=1, ls='--', alpha=0.5)

ytlabels = []
for _, row in sub3.iterrows():
    d = f"{int(row['date_sort'])} BP" if row['date_sort'] > 0 else 'date unknown'
    ytlabels.append(f"{row['Pop2']}  ({d})")
ax.set_yticks(range(len(sub3))); ax.set_yticklabels(ytlabels, fontsize=7.5)
ax.set_xlabel('f4 estimate (±1.96 SE)')
ax.set_title(
    f'f4(Mbuti.DG, Tibetan Pop; Xingyi_EN, YR_MN)\n'
    'Positive = closer to Xingyi_EN  |  Negative = closer to YR_MN\n'
    f'nSNPs ≥ {NSNP_CUTOFF//1000}K; sorted older → younger',
    fontweight='bold'
)
# Legend patches for sub-regions
sp_handles = [
    mpatches.Patch(color='#1A6FBF', label='NE Tibetan (Zongri)'),
    mpatches.Patch(color='#6A3D9A', label='NE Tibetan (Yushu)'),
    mpatches.Patch(color='#33A02C', label='Shannan (Central-S)'),
    mpatches.Patch(color='#E31A1C', label='Chamdo (E Tibet)'),
    mpatches.Patch(color='#FF7F00', label='Shigatse (Central-W)'),
    mpatches.Patch(color='#F4A500', label='Lhasa'),
    mpatches.Patch(color='#B15928', label='Nagqu (Central-N)'),
    mpatches.Patch(color='#1F78B4', label='Ngari (W Tibet)'),
    mpatches.Patch(color='#A6CEE3', label='Nyingchi (SE Tibet)'),
    mpatches.Patch(color='#888888', label='Modern / unlabelled TP'),
]
ax.legend(handles=sp_handles, fontsize=7, loc='lower right', title='Sub-region', title_fontsize=8)
plt.tight_layout()
p = os.path.join(OUT_DIR, f'f4_tibetan_sites_xingyi_en_{DATE_STAMP}.png')
plt.savefig(p, bbox_inches='tight'); plt.close(); print(f"  → {p}")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 4: Tibetan Plateau heatmap (TP Pop2 × selected ref pairs)
# ══════════════════════════════════════════════════════════════════════════════
print("[4/7] Tibetan Plateau heatmap …")

SELECTED_PAIRS = [
    ('Xingyi_EN',       'YR_MN'),
    ('STM_EN',          'Shimao_group'),
    ('Yangshao_UYR',    'Baikal_EN'),
    ('Upper_YR_LN',     'Xingyi_EN'),
    ('YR_MN',           'Baikal_EN'),
    ('Mbc3.5k_G0',      'Shimao_group'),
    ('Mbc3.5k_outlier', 'STM_EN'),
    ('STM_EN',          'Boisman_MN'),
    ('DevilsCave_N',    'Baikal_EN'),
    ('Tianyuan',        'Onge.DG'),
]

# Use pop2_meta dates (includes regex-fallback dates, not just AADR-resolved ones)
tp_date_lookup = dict(zip(meta_df['Pop2'], meta_df['date_bp']))
tp_ancient = sorted(
    [(p, tp_date_lookup.get(p) or 0) for p in tp_pops],
    key=lambda x: -x[1]
)
tp_ancient_pops = [p for p, d in tp_ancient if d > 0]

hm_rows = []
for pop2 in tp_ancient_pops:
    for (p3, p4) in SELECTED_PAIRS:
        row = df[(df['Pop2'] == pop2) & (df['Pop3'] == p3) & (df['Pop4'] == p4)]
        if len(row) == 1:
            hm_rows.append({'Pop2': pop2,
                            'pair': f'{p3}\nvs {p4}',
                            'f4est': row['f4est'].values[0]})
        else:
            row2 = df[(df['Pop2'] == pop2) & (df['Pop3'] == p4) & (df['Pop4'] == p3)]
            if len(row2) == 1:
                hm_rows.append({'Pop2': pop2,
                                'pair': f'{p3}\nvs {p4}',
                                'f4est': -row2['f4est'].values[0]})

hm_df = pd.DataFrame(hm_rows)
if len(hm_df) > 0:
    pivot = hm_df.pivot_table(index='Pop2', columns='pair', values='f4est', aggfunc='mean')
    order = [p for p in tp_ancient_pops if p in pivot.index]
    pivot = pivot.loc[order]

    vmax = max(abs(pivot.values[~np.isnan(pivot.values)]).max(), 0.001)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(max(10, pivot.shape[1]*1.6), max(8, len(pivot)*0.45)))
    sns.heatmap(pivot, ax=ax, cmap='RdBu_r', norm=norm,
                linewidths=0.3, linecolor='white',
                annot=True, fmt='.4f', annot_kws={'size': 6.5},
                cbar_kws={'label': 'f4 estimate', 'shrink': 0.6})
    date_bp_vals = {p: d for p, d in tp_ancient if p in pivot.index and d > 0}
    new_ylabels = [f"{p}  ({int(date_bp_vals.get(p, 0))} BP)" if date_bp_vals.get(p) else f"{p}" for p in pivot.index]
    ax.set_yticklabels(new_ylabels, rotation=0, fontsize=7.5)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha='center', fontsize=8)
    ax.set_title(
        f'f4 Estimates: Tibetan Plateau Ancient Populations\n'
        f'f4(Mbuti.DG, Pop2; Pop3, Pop4) — nSNPs ≥ {NSNP_CUTOFF//1000}K',
        fontweight='bold', pad=10
    )
    ax.set_xlabel('Reference Pop3 vs Pop4 pair', labelpad=8)
    ax.set_ylabel('Pop2 (older → younger)', labelpad=8)
    plt.tight_layout()
    p = os.path.join(OUT_DIR, f'f4_tp_heatmap_{DATE_STAMP}.png')
    plt.savefig(p, bbox_inches='tight'); plt.close(); print(f"  → {p}")
else:
    print("  WARNING: insufficient data for heatmap — skipped.")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 5: Tibetan Plateau vs Himalayan comparison (same reference pairs)
# ══════════════════════════════════════════════════════════════════════════════
print("[5/7] TP vs Himalayan comparison …")

COMPARE_PAIRS = [
    ('Xingyi_EN',    'YR_MN'),
    ('STM_EN',       'Shimao_group'),
    ('Upper_YR_LN',  'Xingyi_EN'),
    ('Mbc3.5k_outlier', 'STM_EN'),
]
# Use ancient TP only (date > 0)
cmp_pops = (
    [(p, 'Tibetan_Plateau') for p in tp_ancient_pops] +
    [(p, 'Himalayan') for p in him_pops if p in df['Pop2'].unique()]
)

cmp_rows = []
for (pop2, grp) in cmp_pops:
    for (p3, p4) in COMPARE_PAIRS:
        row = df[(df['Pop2']==pop2) & (df['Pop3']==p3) & (df['Pop4']==p4)]
        if len(row) == 1:
            cmp_rows.append({'Pop2': pop2, 'Group': grp,
                             'pair': f'{p3}\nvs {p4}',
                             'f4est': row['f4est'].values[0],
                             'se': row['se'].values[0],
                             'z': row['z'].values[0]})

cmp_df = pd.DataFrame(cmp_rows)
if len(cmp_df) > 0:
    pairs_avail = cmp_df['pair'].unique()
    ncols = len(pairs_avail)
    fig, axes = plt.subplots(1, ncols, figsize=(5.5*ncols, max(8, len(cmp_pops)*0.3)),
                              sharey=False)
    if ncols == 1: axes = [axes]
    grp_colors = {'Tibetan_Plateau': '#2166AC', 'Himalayan': '#74ADD1'}

    for ax, pair in zip(axes, pairs_avail):
        sub = cmp_df[cmp_df['pair'] == pair].sort_values(['Group', 'f4est'])
        for i, (_, row) in enumerate(sub.iterrows()):
            c = grp_colors.get(row['Group'], '#888')
            ax.errorbar(row['f4est'], i, xerr=1.96*row['se'],
                        fmt='o', color=c, ecolor=c, ms=5, elinewidth=1.2, capsize=2.5)
        ax.axvline(0, color='black', lw=0.7)
        ax.set_yticks(range(len(sub))); ax.set_yticklabels(sub['Pop2'], fontsize=7)
        ax.set_title(pair.replace('\n', ' '), fontsize=9, fontweight='bold')
        ax.set_xlabel('f4 estimate (±1.96 SE)', fontsize=8)

    patches = [mpatches.Patch(color=c, label=g) for g, c in grp_colors.items()]
    fig.legend(handles=patches, loc='upper right', fontsize=9)
    fig.suptitle(
        f'f4(Mbuti.DG, Pop2; Pop3, Pop4)\nTibetan Plateau vs Himalayan/Mustang — nSNPs ≥ {NSNP_CUTOFF//1000}K',
        fontsize=12, fontweight='bold'
    )
    plt.tight_layout()
    p = os.path.join(OUT_DIR, f'f4_tibet_mustang_comparison_{DATE_STAMP}.png')
    plt.savefig(p, bbox_inches='tight'); plt.close(); print(f"  → {p}")
else:
    print("  WARNING: no data for TP vs Himalayan comparison — skipped.")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 6: Temporal gradient — dated TP pops across three reference pairs
# ══════════════════════════════════════════════════════════════════════════════
print("[6/7] Temporal gradient …")

TEMP_PAIRS = [
    ('Xingyi_EN',       'YR_MN'),
    ('Upper_YR_LN',     'Xingyi_EN'),
    ('STM_EN',          'Shimao_group'),
]
# All TP pops with known date (AADR or regex fallback), sorted oldest→newest
dated_tp = sorted(
    [(p, tp_date_lookup[p]) for p in tp_pops
     if tp_date_lookup.get(p) and tp_date_lookup[p] > 0],
    key=lambda x: -x[1]
)
# Colour ramp: use a continuous blue-based palette by date
cmap_temp = plt.cm.get_cmap('plasma_r', len(dated_tp))
pop_color_map = {pop: cmap_temp(i / max(len(dated_tp)-1, 1))
                 for i, (pop, _) in enumerate(dated_tp)}

fig, axes = plt.subplots(1, len(TEMP_PAIRS), figsize=(5.5*len(TEMP_PAIRS), max(9, len(dated_tp)*0.38)))

for ax, (p3, p4) in zip(axes, TEMP_PAIRS):
    rows = []
    for (pop2, date_bp) in dated_tp:
        row = df[(df['Pop2']==pop2) & (df['Pop3']==p3) & (df['Pop4']==p4)]
        if len(row) == 1:
            rows.append({'Pop2': pop2, 'date_bp': date_bp,
                         'f4est': row['f4est'].values[0],
                         'se':    row['se'].values[0],
                         'z':     row['z'].values[0]})
    if not rows: continue
    sub = pd.DataFrame(rows).sort_values('date_bp', ascending=False)

    for i, (_, r) in enumerate(sub.iterrows()):
        c = pop_color_map.get(r['Pop2'], '#888')
        sig = abs(r['z']) >= 3
        ax.errorbar(r['f4est'], i, xerr=1.96*r['se'],
                    fmt='*' if sig else 'o', color=c, ecolor=c,
                    ms=8 if sig else 5, elinewidth=1.4, capsize=2.5, zorder=3)

    ax.axvline(0, color='black', lw=0.8)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels([f"{r['Pop2']} ({int(r['date_bp'])} BP)" for _, r in sub.iterrows()],
                        fontsize=7.5)
    ax.set_xlabel('f4 estimate (±1.96 SE)')
    ax.set_title(f'Pop3: {p3}\nPop4: {p4}', fontsize=9, fontweight='bold')
    ax.text(0.02, 0.01, '★ = |Z| ≥ 3', transform=ax.transAxes, fontsize=8, color='black')

fig.suptitle(
    f'Temporal Gradient: f4(Mbuti.DG, Tibetan Pop; Pop3, Pop4)\n'
    f'Older → younger (top → bottom) | nSNPs ≥ {NSNP_CUTOFF//1000}K',
    fontsize=12, fontweight='bold'
)
plt.tight_layout()
p = os.path.join(OUT_DIR, f'f4_temporal_gradient_{DATE_STAMP}.png')
plt.savefig(p, bbox_inches='tight'); plt.close(); print(f"  → {p}")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 7: Mabuco / Sichuan affinity forest plot
# ══════════════════════════════════════════════════════════════════════════════
print("[7/7] Mabuco / Sichuan affinity …")

avail_mbc = [p for p in mbc_pops if p in df['Pop2'].unique()]
MBC_PAIRS = [
    ('Upper_YR_LN',     'Xingyi_EN'),
    ('Xingyi_EN',       'YR_MN'),
    ('STM_EN',          'Shimao_group'),
    ('Mbc3.5k_G0',      'Shimao_group'),
    ('Mbc3.5k_outlier', 'STM_EN'),
]

mbc_rows = []
for pop2 in avail_mbc:
    for (p3, p4) in MBC_PAIRS:
        row = df[(df['Pop2']==pop2) & (df['Pop3']==p3) & (df['Pop4']==p4)]
        if len(row) == 1:
            d_bp = pid_date_map.get(pop2, None)
            mbc_rows.append({'Pop2': pop2,
                             'pair': f'{p3} vs {p4}',
                             'date_bp': d_bp,
                             'f4est': row['f4est'].values[0],
                             'se': row['se'].values[0],
                             'z': row['z'].values[0]})

mbc_df = pd.DataFrame(mbc_rows)
if len(mbc_df) > 0:
    pairs_avail = mbc_df['pair'].unique()
    fig, axes = plt.subplots(
        1, len(pairs_avail),
        figsize=(4.8*len(pairs_avail), max(6, len(avail_mbc)*0.45)),
        sharey=True
    )
    if len(pairs_avail) == 1: axes = [axes]

    # Colour by pop prefix
    def mbc_color(pop):
        if pop.startswith('MabucoE'):       return '#D73027'
        if pop.startswith('Mbc3.5k'):       return '#FC8D59'
        if pop.startswith('Mbc4.4k'):       return '#74ADD1'
        if pop.startswith('Mbc4k'):         return '#4575B4'
        if pop.startswith('Xingyi_EN'):     return '#006837'
        if pop.startswith('Xingyi_BA'):     return '#31A354'
        if pop.startswith('Xingyi_LN'):     return '#78C679'
        if pop.startswith('DHL'):           return '#A6761D'
        if pop.startswith('Haneyi'):        return '#984EA3'
        if pop.startswith('BGD'):           return '#762A83'
        return '#888'

    first_ax = axes[0]
    for ax, pair in zip(axes, pairs_avail):
        sub = mbc_df[mbc_df['pair'] == pair].copy()
        # Sort by date if available, else by f4est
        sub['ds'] = sub['date_bp'].fillna(-1)
        sub = sub.sort_values('ds', ascending=False)
        for i, (_, row) in enumerate(sub.iterrows()):
            c = mbc_color(row['Pop2'])
            sig = abs(row['z']) >= 3
            ax.errorbar(row['f4est'], i, xerr=1.96*row['se'],
                        fmt='*' if sig else 'o', color=c, ecolor=c,
                        ms=8 if sig else 6, elinewidth=1.3, capsize=2.5)
        ax.axvline(0, color='black', lw=0.8)
        ax.set_yticks(range(len(sub)))
        ylabels = []
        for _, r in sub.iterrows():
            d = f"  ({int(r['ds'])} BP)" if r['ds'] > 0 else ''
            ylabels.append(r['Pop2'] + d)
        ax.set_yticklabels(ylabels, fontsize=7.5)
        ax.set_title(pair, fontsize=9, fontweight='bold')
        ax.set_xlabel('f4 estimate (±1.96 SE)', fontsize=8)
        ax.text(0.02, 0.01, '★ = |Z| ≥ 3', transform=ax.transAxes, fontsize=8)

    fig.suptitle(
        f'f4(Mbuti.DG, Sichuan/Mabuco Pop; Pop3, Pop4)\nnSNPs ≥ {NSNP_CUTOFF//1000}K',
        fontsize=12, fontweight='bold'
    )
    plt.tight_layout()
    p = os.path.join(OUT_DIR, f'f4_mbc_sichuan_affinity_{DATE_STAMP}.png')
    plt.savefig(p, bbox_inches='tight'); plt.close(); print(f"  → {p}")
else:
    print("  WARNING: no Mabuco/Sichuan data found — skipped.")


# ══════════════════════════════════════════════════════════════════════════════
# Final report
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*65)
print("FINAL REPORT")
print("="*65)
generated = sorted(f for f in os.listdir(OUT_DIR) if f.endswith('.png'))
print(f"\nGenerated figures ({len(generated)}) in {OUT_DIR}:")
for f in generated:
    print(f"  {f}")

print(f"\nData summary:")
print(f"  f4 raw rows:             {n_raw:,}")
print(f"  After nSNPs ≥ {NSNP_CUTOFF//1000}K:     {len(df):,}  (removed {n_raw-len(df):,})")
print(f"  Z-score NA rows dropped: {len(df) - len(df_noNA)}")

print(f"\nPop2 groups matched to date via .ind → AADR:")
print(f"  With date resolved:      {matched_count}")
print(f"  Date unknown (metadata_not_found): {unmatched_count}")

print(f"\nRegional breakdown of Pop2 in data:")
for reg, grp in meta_df.groupby('region'):
    n_dated = grp['date_bp'].notna().sum()
    print(f"  {reg:<20} {len(grp):>4} pops  ({n_dated} with date)")

print(f"\nMetadata files used:")
print(f"  {IND_PATH}  — IID→PID mapping")
print(f"  {AADR_PATH} — date lookup via Genetic ID")
print(f"  {SI_PATH}   — period/cal labels (supplementary)")
print(f"\nMatching strategy: IID from .ind → AADR 'Genetic ID' column → date mean in BP averaged per PID")
print(f"\nCaveats:")
print(f"  • Groups without names ending in a date suffix are NOT assumed to be modern;")
print(f"    their dates are resolved via AADR lookup.")
print(f"  • 'Outlier' sub-groups (e.g. Zongri4.5k_o1) may represent excluded individuals.")
print(f"  • Date is mean across all AADR-matched individuals in the group.")
print(f"  • Some populations (SDLG, GBSL, QLSZ, PLDW etc.) returned no date from AADR;")
print(f"    they appear in heatmaps/comparison plots but are excluded from temporal-gradient figures.")
print("\nDone.")
