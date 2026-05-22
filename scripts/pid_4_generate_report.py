import os
import pandas as pd
import markdown
import base64
import glob

BASE_DIR = "/home/duri_bae/popgen_lab"
OUT_DIR = os.path.join(BASE_DIR, "output/adaptive_sweep/PID_comparison")
FIG_DIR = os.path.join(OUT_DIR, "figures")

MD_FILE = os.path.join(OUT_DIR, "new_pid_adaptive_sweep_comparison_report_260522.md")
HTML_FILE = os.path.join(OUT_DIR, "new_pid_adaptive_sweep_comparison_report_260522.html")

presence = pd.read_csv(os.path.join(OUT_DIR, "new_pid_presence_summary_260522.csv"))
candidates = pd.read_csv(os.path.join(OUT_DIR, "new_pid_characteristic_snp_candidates_260522.csv"))

md = []
md.append("# Adaptive Sweep: New PID Comparison Report")
md.append("**Date**: 2026-05-22")
md.append("**Task Objective**: Extend the previous adaptive sweep analysis by treating a new set of 18 PIDs as independent populations. We compare these new PIDs against the previously established early Tibetan targets to identify PID-specific fixed, high-frequency, and highly differentiated SNPs.")

md.append("## 1. Input Files & Data Reuse")
md.append("- **Input Genotypes**: `/raw/Tibetan.fstat.1240K.260508.bed/bim/fam`")
md.append("- **Previous Targets Inferred**: `Zongri4.7k, Mbc3.5k_G0, Mbc4.4k_G2, Zongri5.1k, GBSL, GBSL_old1, Shannan3k, Mbc4k_G2, SUI001.WGS, Yushu2.8k, Zongri4.1k, Chamdo2.8k_1, Mbc3.5k_outlier, Mbc4.4k_G0, Zongri4.5k, GBSL_old2, Mbc4.4k_G1`")
md.append("- **Inference Method**: The previous target populations and validated known SNPs were dynamically extracted from the text and tables of the previous `/output/adaptive_sweep/` reports.")

md.append("## 2. PID Presence Summary")
md.append(presence.to_markdown(index=False))

md.append("## 3. Methods")
md.append("- **Genotype Extraction**: Pure Python binary matrix extraction from the `.bed` file to bypass missing dependencies.")
md.append("- **IID-based Region/Period Assignment**: Matched exactly against `SampleInfo_260522.csv` and `v66.1240K.aadr.PUB.csv` using the internal individual IDs (IID) rather than guessing population labels.")
md.append("- **SNP-Gene Mapping**: Strict reuse of the exact `snp_label` and mapping rationale established in the previous adaptive sweep analysis.")
md.append("- **Contrast Calculation**: Pooled Target AF (aggregated across the 17 early Tibetan populations) was subtracted from each new PID's independent AF.")

md.append("## 4. PID-Specific Characteristic SNP Candidates")
md.append(f"We identified **{len(candidates)}** characteristic SNP occurrences across the 18 PIDs that exhibit extreme divergence or fixation. Below is the top 15 most differentiated list:")
md.append(candidates.sort_values(by="AF_Difference", key=abs, ascending=False).head(15).to_markdown(index=False))

md.append("## 5. Visualizations")
# Gather all pngs
pngs = glob.glob(os.path.join(FIG_DIR, "*.png"))
pngs.sort()

md_for_md = md.copy()
md_for_html = md.copy()

for p in pngs:
    fname = os.path.basename(p)
    title = fname.replace("_260522.png", "").replace("_", " ").title()
    
    md_for_md.append(f"### {title}")
    md_for_md.append(f"![{title}](figures/{fname})")
    
    with open(p, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode()
    
    md_for_html.append(f"### {title}")
    md_for_html.append(f"<img alt=\"{title}\" src=\"data:image/png;base64,{encoded_string}\" />")

md_for_md.append("## 6. Limitations & Caveats")
md.append("1. **Missingness & Ascertainment**: Ancient DNA suffers heavily from pseudo-haploid sampling and low-coverage dropout. Allele frequencies based on n < 5 should be strictly interpreted as `caution_low_sample_size`.")
md.append("2. **Selection Inference**: High allele frequency alone is **not** definitive proof of positive selection. Genetic drift, population bottlenecks, and capture bias can equally cause fixation.")
md.append("3. **Independent Grouping**: The PIDs were analyzed entirely independently. We did not pool them into a meta-reference group.")

md_for_md.append("## 7. Output Manifest")
md_for_md.append("- `/new_pid_presence_summary_260522.csv`\n- `/new_pid_known_snp_frequencies_260522.csv`\n- `/previous_targets_vs_each_new_pid_snp_contrast_260522.csv`\n- `/new_pid_characteristic_snp_candidates_260522.csv`")
md_for_html.append("## 6. Limitations & Caveats")
md_for_html.append("1. **Missingness & Ascertainment**: Ancient DNA suffers heavily from pseudo-haploid sampling and low-coverage dropout. Allele frequencies based on n < 5 should be strictly interpreted as `caution_low_sample_size`.")
md_for_html.append("2. **Selection Inference**: High allele frequency alone is **not** definitive proof of positive selection. Genetic drift, population bottlenecks, and capture bias can equally cause fixation.")
md_for_html.append("3. **Independent Grouping**: The PIDs were analyzed entirely independently. We did not pool them into a meta-reference group.")
md_for_html.append("## 7. Output Manifest")
md_for_html.append("- `/new_pid_presence_summary_260522.csv`\n- `/new_pid_known_snp_frequencies_260522.csv`\n- `/previous_targets_vs_each_new_pid_snp_contrast_260522.csv`\n- `/new_pid_characteristic_snp_candidates_260522.csv`")

md_text_for_md = "\n\n".join(md_for_md)
with open(MD_FILE, "w", encoding="utf-8") as f:
    f.write(md_text_for_md)

html_body = markdown.markdown("\n\n".join(md_for_html), extensions=['tables'])

html_full = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>PID Comparison Sweep Report</title>
    <style>
        :root {{
            --primary-color: #2c3e50;
            --bg-color: #f8f9fa;
        }}
        body {{ 
            font-family: 'Inter', 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
            line-height: 1.6; 
            margin: 0; 
            padding: 40px 20px; 
            background-color: var(--bg-color);
            color: #333;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: #fff;
            padding: 40px;
            border-radius: 8px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.05);
        }}
        h1, h2, h3 {{ color: var(--primary-color); border-bottom: 2px solid #eee; padding-bottom: 10px; margin-top: 40px; }}
        h1 {{ margin-top: 0; }}
        .table-wrapper {{
            overflow-x: auto;
            margin-bottom: 20px;
            border-radius: 4px;
            border: 1px solid #e0e0e0;
        }}
        table {{ 
            border-collapse: collapse; 
            width: 100%; 
            white-space: nowrap; 
        }}
        th, td {{ border: 1px solid #e0e0e0; padding: 12px 15px; text-align: left; }}
        th {{ background-color: #f1f3f5; color: #495057; font-weight: 600; position: sticky; top: 0; }}
        tr:nth-child(even) {{ background-color: #fafafa; }}
        img {{ max-width: 100%; height: auto; border-radius: 6px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin: 20px 0; display: block; }}
    </style>
</head>
<body>
<div class="container">
{html_body}
</div>
<script>
    document.querySelectorAll('table').forEach(function(table) {{
        var wrapper = document.createElement('div');
        wrapper.className = 'table-wrapper';
        table.parentNode.insertBefore(wrapper, table);
        wrapper.appendChild(table);
    }});
</script>
</body>
</html>
"""

with open(HTML_FILE, "w", encoding="utf-8") as f:
    f.write(html_full)

print(f"Report generated successfully at {HTML_FILE}")
