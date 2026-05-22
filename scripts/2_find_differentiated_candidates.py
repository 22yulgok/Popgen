import os
import csv
import numpy as np
import collections

OUTPUT_DIR = "/home/duri_bae/popgen_lab/output/adaptive_sweep"
BED_FILE = "/home/duri_bae/popgen_lab/raw/Tibetan.fstat.1240K.260508.bed"
BIM_FILE = "/home/duri_bae/popgen_lab/raw/Tibetan.fstat.1240K.260508.bim"
FAM_FILE = "/home/duri_bae/popgen_lab/raw/Tibetan.fstat.1240K.260508.fam"
OUT_CSV = os.path.join(OUTPUT_DIR, "high_frequency_polymorphic_candidates_260522.csv")

TARGET_POPS = {'GBSL', 'GBSL_old1', 'GBSL_old2', 'Yushu2.8k', 'Zongri5.1k', 'Zongri4.7k', 'Zongri4.1k', 'Zongri4.5k', 'Shannan3k', 'Chamdo2.8k_1', 'Mbc3.5k_G0', 'Mbc3.5k_outlier', 'Mbc4.4k_G0', 'Mbc4.4k_G1', 'Mbc4.4k_G2', 'Mbc4k_G2', 'SUI001.WGS'}

print("Reading FAM...")
samples = []
target_indices = []
comp_indices = []
pop_to_indices = collections.defaultdict(list)

with open(FAM_FILE, "r") as f:
    for i, line in enumerate(f):
        parts = line.strip().split()
        if len(parts) >= 2:
            fid, iid = parts[0], parts[1]
            samples.append({"FID": fid, "IID": iid})
            pop_to_indices[fid].append(i)
            if fid in TARGET_POPS:
                target_indices.append(i)
            else:
                comp_indices.append(i)

N = len(samples)
target_indices = np.array(target_indices)
comp_indices = np.array(comp_indices)

target_pops_filtered = {p: np.array(pop_to_indices[p]) for p in TARGET_POPS if len(pop_to_indices[p]) >= 3}

print("Reading BIM...")
bim_snps = []
with open(BIM_FILE, "r") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 6:
            bim_snps.append({"chr": parts[0], "rsid": parts[1], "pos": parts[3], "a1": parts[4], "a2": parts[5]})

num_snps = len(bim_snps)
bytes_per_snp = int(np.ceil(N / 4.0))

print("Reading BED chunks...")
candidates = []

with open(BED_FILE, "rb") as f:
    magic = f.read(3)
    if magic != b"\x6c\x1b\x01":
        raise ValueError("Invalid BED magic")
    bed_data = np.fromfile(f, dtype=np.uint8)

# Trim if there's extra bytes (some PLINK files have padding, but usually exact)
expected_bytes = num_snps * bytes_per_snp
if len(bed_data) > expected_bytes:
    bed_data = bed_data[:expected_bytes]
elif len(bed_data) < expected_bytes:
    print(f"Warning: expected {expected_bytes} bytes but got {len(bed_data)}")

bed_matrix = bed_data.reshape((num_snps, bytes_per_snp))

CHUNK_SIZE = 50000
shifts = np.array([0, 2, 4, 6], dtype=np.uint8)

for start in range(0, num_snps, CHUNK_SIZE):
    end = min(start + CHUNK_SIZE, num_snps)
    chunk = bed_matrix[start:end]
    
    # Fast bit unpacking
    expanded = (chunk[:, :, None] >> shifts) & 3
    genos = expanded.reshape(chunk.shape[0], -1)[:, :N]
    
    valid_mask = (genos != 1)
    a2_copies = np.zeros_like(genos)
    a2_copies[genos == 2] = 1
    a2_copies[genos == 3] = 2
    
    tgt_valid = valid_mask[:, target_indices]
    tgt_a2 = a2_copies[:, target_indices]
    tgt_called = tgt_valid.sum(axis=1)
    tgt_a2_sum = (tgt_a2 * tgt_valid).sum(axis=1)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        tgt_af = np.where(tgt_called > 0, tgt_a2_sum / (2.0 * tgt_called), np.nan)
        
    cmp_valid = valid_mask[:, comp_indices]
    cmp_a2 = a2_copies[:, comp_indices]
    cmp_called = cmp_valid.sum(axis=1)
    cmp_a2_sum = (cmp_a2 * cmp_valid).sum(axis=1)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        cmp_af = np.where(cmp_called > 0, cmp_a2_sum / (2.0 * cmp_called), np.nan)
        
    for j in range(end - start):
        if tgt_called[j] < 5 or cmp_called[j] < 20: continue
        
        t_af = tgt_af[j]
        c_af = cmp_af[j]
        
        if np.isnan(t_af) or np.isnan(c_af): continue
        
        diff = t_af - c_af
        is_candidate = False
        cat = ""
        
        if abs(diff) >= 0.7:
            is_candidate = True
            cat = "differentiated_candidate"
        elif t_af >= 0.95 and c_af <= 0.5:
            is_candidate = True
            cat = "high_frequency_candidate"
        elif t_af <= 0.05 and c_af >= 0.5:
            is_candidate = True
            cat = "high_frequency_candidate"
            
        if is_candidate:
            idx = start + j
            b = bim_snps[idx]
            
            # Check fixation in specific target populations
            g_row = genos[j]
            v_row = valid_mask[j]
            a_row = a2_copies[j]
            
            fixed_pops = []
            for p, p_idx in target_pops_filtered.items():
                p_v = v_row[p_idx]
                p_a = a_row[p_idx]
                p_call = p_v.sum()
                if p_call >= 3:
                    p_af = (p_a * p_v).sum() / (2.0 * p_call)
                    if t_af > 0.5 and p_af >= 0.95:
                        fixed_pops.append(p)
                    elif t_af < 0.5 and p_af <= 0.05:
                        fixed_pops.append(p)
            
            if len(fixed_pops) > 0 and abs(diff) >= 0.5:
                cat = "fixed_candidate"
                
            candidates.append({
                "rsID": b["rsid"],
                "chr": b["chr"],
                "pos": b["pos"],
                "A1": b["a1"],
                "A2": b["a2"],
                "Target_AF_A2": round(t_af, 4),
                "Comp_AF_A2": round(c_af, 4),
                "Category": cat,
                "Fixed_Target_Pops": ";".join(fixed_pops),
                "Annotation": "No verified selection-related evidence found for this SNP or nearby gene."
            })

# Keep top 300 candidates
candidates.sort(key=lambda x: abs(x["Target_AF_A2"] - x["Comp_AF_A2"]), reverse=True)
candidates = candidates[:300]

with open(OUT_CSV, "w") as f:
    writer = csv.DictWriter(f, fieldnames=["rsID", "chr", "pos", "A1", "A2", "Target_AF_A2", "Comp_AF_A2", "Category", "Fixed_Target_Pops", "Annotation"])
    writer.writeheader()
    for c in candidates:
        writer.writerow(c)

print(f"Saved {len(candidates)} candidate SNPs.")
