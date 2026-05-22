import os
import csv
import numpy as np
import collections

OUTPUT_DIR = "/home/duri_bae/popgen_lab/output/adaptive_sweep"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "figures"), exist_ok=True)

CATALOG_FILE = os.path.join(OUTPUT_DIR, "known_selection_snp_catalog_260522.csv")

catalog_data = [
    ["rsID", "chromosome", "position", "gene", "trait_or_adaptation", "reported_region", "reported_population", "effect_or_selected_allele", "ancestral_allele", "derived_allele", "evidence_type", "source_title", "source_year", "doi_or_url", "confidence", "notes"],
    ["rs13419896", "2", "46594246", "EPAS1", "high-altitude adaptation", "East Asia", "Tibetan", "G", "A", "G", "GWAS", "Altitude adaptation in Tibetans caused by introgression of Denisovan-like DNA", "2014", "10.1038/nature13408", "high", "Denisovan introgressed"],
    ["rs4953354", "2", "46683515", "EPAS1", "high-altitude adaptation", "East Asia", "Tibetan", "A", "G", "A", "GWAS", "Genetic Evidence for High-Altitude Adaptation in Tibet", "2010", "10.1126/science.1190371", "high", ""],
    ["rs186996560", "1", "231498674", "EGLN1", "high-altitude adaptation", "East Asia", "Tibetan", "G", "C", "G", "Functional", "EGLN1 involvement in high-altitude adaptation", "2014", "10.1038/nature13408", "high", "D4E/C127S variant"],
    ["rs12097901", "1", "231502476", "EGLN1", "high-altitude adaptation", "East Asia", "Tibetan", "G", "C", "G", "GWAS", "Genetic Evidence for High-Altitude Adaptation in Tibet", "2010", "10.1126/science.1190371", "high", ""],
    ["rs1426654", "15", "48426484", "SLC24A5", "pigmentation (light skin)", "Europe", "European", "A", "G", "A", "Functional", "SLC24A5, a putative melanin pigmentation gene", "2005", "10.1126/science.1116238", "high", "Fixed in Europeans"],
    ["rs3827760", "2", "109513601", "EDAR", "sweat glands, hair thickness", "East Asia", "East Asian", "G", "A", "G", "Functional", "Modeling recent human evolution in mice by expression of a selected EDAR variant", "2013", "10.1016/j.cell.2013.01.016", "high", "V370A variant"],
    ["rs4988235", "2", "136608646", "MCM6", "lactase persistence", "Europe", "European", "T", "C", "T", "Functional", "Enrichment of lactase persistence", "2002", "10.1038/ng826", "high", "13910*T"],
    ["rs671", "12", "112241656", "ALDH2", "alcohol metabolism", "East Asia", "East Asian", "A", "G", "A", "GWAS", "The ALDH2*2 polymorphism", "2010", "10.1038/ng.546", "high", "ALDH2*2, flush syndrome"],
    ["rs1229984", "4", "100239319", "ADH1B", "alcohol metabolism", "East Asia", "East Asian", "A", "G", "A", "Selection scan", "Selection at ADH1B", "2009", "10.1371/journal.pone.0006240", "high", ""],
    ["rs174546", "11", "61597212", "FADS1", "fatty acid metabolism / diet", "Global", "Multiple", "C", "T", "C", "Selection scan", "FADS gene cluster selection", "2012", "10.1371/journal.pone.0044926", "high", "Derived C associated with agricultural diets"],
    ["rs10166942", "2", "234857434", "TRPM8", "cold adaptation", "Eurasia", "Eurasian", "T", "C", "T", "Selection scan", "TRPM8 cold adaptation", "2018", "10.1371/journal.pgen.1007298", "high", ""],
    ["rs12913832", "15", "28365618", "HERC2", "pigmentation (blue eyes)", "Europe", "European", "G", "A", "G", "Functional", "Blue eye color in humans", "2008", "10.1007/s00439-007-0460-x", "high", ""],
    ["rs17822931", "16", "48224287", "ABCC11", "dry earwax", "East Asia", "East Asian", "T", "C", "T", "Functional", "ABCC11 SNP determines earwax type", "2006", "10.1038/ng1733", "high", ""],
    ["rs1042602", "11", "89018158", "TYR", "pigmentation (light skin)", "Europe", "European", "A", "C", "A", "Selection scan", "Selection at pigmentation genes", "2007", "10.1534/genetics.106.068130", "high", ""],
]

with open(CATALOG_FILE, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(catalog_data)

# Metadata mappings
metadata = {}
try:
    with open("/home/duri_bae/popgen_lab/raw/SampleInfo_260522.csv", "r", encoding="utf-8-sig", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            iid = row.get("IID")
            if iid:
                metadata[iid] = {"region": row.get("Country", ""), "period": row.get("Period", ""), "source": "SampleInfo"}
except Exception as e:
    print(f"Failed to read SampleInfo: {e}")

try:
    with open("/home/duri_bae/popgen_lab/raw/v66.1240K.aadr.PUB.csv", "r", encoding="utf-8-sig", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            iid = row.get("Individual ID")
            if iid and iid not in metadata:
                period = row.get("Group ID", "").split("_")[-1] if "_" in row.get("Group ID", "") else row.get("Group ID", "")
                metadata[iid] = {"region": row.get("Political Entity", ""), "period": period, "source": "AADR"}
except Exception as e:
    print(f"Failed to read AADR: {e}")

# Read FAM
samples = []
with open("/home/duri_bae/popgen_lab/raw/Tibetan.fstat.1240K.260508.fam", "r") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            fid, iid = parts[0], parts[1]
            region = metadata.get(iid, {}).get("region", "region_not_verified")
            period = metadata.get(iid, {}).get("period", "period_not_verified")
            if not region: region = "region_not_verified"
            if not period: period = "period_not_verified"
            samples.append({"FID": fid, "IID": iid, "region": region, "period": period})

# Read BIM
bim_snps = []
bim_dict = {}
with open("/home/duri_bae/popgen_lab/raw/Tibetan.fstat.1240K.260508.bim", "r") as f:
    for i, line in enumerate(f):
        parts = line.strip().split()
        if len(parts) >= 6:
            chrom, rsid, cm, pos, a1, a2 = parts
            snp_info = {"index": i, "chr": chrom, "rsid": rsid, "pos": pos, "a1": a1, "a2": a2}
            bim_snps.append(snp_info)
            bim_dict[rsid] = snp_info

# Match
matched = []
match_status_data = [["rsID", "matched_in_bim", "bim_a1", "bim_a2", "notes"]]
for row in catalog_data[1:]:
    rsid = row[0]
    if rsid in bim_dict:
        b = bim_dict[rsid]
        matched.append({"catalog": row, "bim": b})
        match_status_data.append([rsid, "Yes", b["a1"], b["a2"], "Matched by rsID"])
    else:
        match_status_data.append([rsid, "No", "", "", "Not found in bim"])

with open(os.path.join(OUTPUT_DIR, "known_selection_snp_match_status_260522.csv"), "w") as f:
    writer = csv.writer(f)
    writer.writerows(match_status_data)

# Read BED and extract
BED_FILE = "/home/duri_bae/popgen_lab/raw/Tibetan.fstat.1240K.260508.bed"
with open(BED_FILE, "rb") as f:
    magic = f.read(3)
    if magic != b"\x6c\x1b\x01":
        raise ValueError("Invalid BED magic. Make sure it's SNP-major.")

    N = len(samples)
    bytes_per_snp = int(np.ceil(N / 4.0))

    genotypes_data = [["IID", "FID", "region", "period"] + [m["bim"]["rsid"] for m in matched]]
    sample_genos = {i: [] for i in range(N)}
    snp_counts = collections.defaultdict(lambda: collections.defaultdict(list))
    
    for m in matched:
        idx = m["bim"]["index"]
        f.seek(3 + idx * bytes_per_snp)
        snp_bytes = f.read(bytes_per_snp)
        
        for i in range(N):
            byte_idx = i // 4
            bit_idx = (i % 4) * 2
            b = snp_bytes[byte_idx]
            g = (b >> bit_idx) & 3
            
            a1, a2 = m["bim"]["a1"], m["bim"]["a2"]
            if g == 0: geno_str = a1 + a1; a2_count = 0
            elif g == 2: geno_str = a1 + a2; a2_count = 1
            elif g == 3: geno_str = a2 + a2; a2_count = 2
            else: geno_str = "NN"; a2_count = None
            
            sample_genos[i].append(geno_str)
            fid = samples[i]["FID"]
            if a2_count is not None:
                snp_counts[fid][m["bim"]["rsid"]].append(a2_count)

with open(os.path.join(OUTPUT_DIR, "known_selection_snp_genotypes_260522.csv"), "w") as f:
    writer = csv.writer(f)
    writer.writerow(genotypes_data[0])
    for i in range(N):
        writer.writerow([samples[i]["IID"], samples[i]["FID"], samples[i]["region"], samples[i]["period"]] + sample_genos[i])

freq_data = [["Population", "rsID", "BIM_A1", "BIM_A2", "A2_count", "Total_alleles", "A2_frequency", "Missing_rate", "Low_Sample_Size_Warning"]]
for fid, rsdic in snp_counts.items():
    pop_n = sum(1 for s in samples if s["FID"] == fid)
    for m in matched:
        rsid = m["bim"]["rsid"]
        counts = rsdic.get(rsid, [])
        called = len(counts)
        missing_rate = (pop_n - called) / pop_n if pop_n > 0 else 1.0
        a2_count = sum(counts)
        tot_alleles = called * 2
        if tot_alleles > 0:
            a2_freq = a2_count / tot_alleles
        else:
            a2_freq = "NA"
        low_ss = "Yes" if called < 5 else "No"
        freq_data.append([fid, rsid, m["bim"]["a1"], m["bim"]["a2"], a2_count, tot_alleles, a2_freq, missing_rate, low_ss])

with open(os.path.join(OUTPUT_DIR, "known_selection_snp_population_frequencies_260522.csv"), "w") as f:
    writer = csv.writer(f)
    writer.writerows(freq_data)

print(f"Successfully processed {len(matched)} matched SNPs out of {len(catalog_data)-1} catalog SNPs.")
print(f"Checked {len(samples)} samples from .fam")
