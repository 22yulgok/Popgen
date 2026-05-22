import os
import pandas as pd
import markdown
import base64
import glob

OUTPUT_DIR = "/home/duri_bae/popgen_lab/output/adaptive_sweep"
FIG_DIR = os.path.join(OUTPUT_DIR, "figures")
MD_FILE = os.path.join(OUTPUT_DIR, "adaptive_sweep_report_260522.md")
HTML_FILE = os.path.join(OUTPUT_DIR, "adaptive_sweep_report_260522.html")

TARGET_POPS = {'GBSL', 'GBSL_old1', 'GBSL_old2', 'Yushu2.8k', 'Zongri5.1k', 'Zongri4.7k', 'Zongri4.1k', 'Zongri4.5k', 'Shannan3k', 'Chamdo2.8k_1', 'Mbc3.5k_G0', 'Mbc3.5k_outlier', 'Mbc4.4k_G0', 'Mbc4.4k_G1', 'Mbc4.4k_G2', 'Mbc4k_G2', 'SUI001.WGS'}

match_status = pd.read_csv(os.path.join(OUTPUT_DIR, "known_selection_snp_match_status_260522.csv"))
freqs = pd.read_csv(os.path.join(OUTPUT_DIR, "known_selection_snp_population_frequencies_260522.csv"))
candidates = pd.read_csv(os.path.join(OUTPUT_DIR, "high_frequency_polymorphic_candidates_260522.csv"))

md = []
md.append("# Adaptive Sweep and Selection Analysis Report")
md.append("**Date**: 2026-05-22  <br>**Task Objective**: Investigate adaptive sweep and selection-related SNP patterns in the target early Tibetan Plateau groups and surrounding reference populations.")

md.append("## 1. Input Files & Populations")
md.append("- **Input Files**: `/raw/Tibetan.fstat.1240K.260508.bed/bim/fam`\n- **Target Populations Found**: " + ", ".join(TARGET_POPS))

md.append("## 2. Known SNP Catalog Summary")
md.append(match_status.to_markdown(index=False))

md.append("### Key Genes Investigated")
md.append("- **EPAS1 / EGLN1**: Key hypoxia pathway genes strongly associated with high-altitude adaptation in Tibetan populations.")
md.append("- **SLC24A5 / OCA2 / TYR**: Genes driving lighter skin and eye pigmentation, typically selected in higher latitudes (e.g., Europeans).")
md.append("- **EDAR**: The V370A variant is associated with increased sweat gland density and thicker hair in East Asians.")
md.append("- **LCT (MCM6)**: The 13910*T variant is the primary driver of lactase persistence in European populations.")
md.append("- **ALDH2 / ADH1B**: Variants conferring the 'alcohol flush' reaction, heavily selected in East Asian populations.")
md.append("- **TRPM8**: A cold-sensing receptor variant implicated in cold adaptation and migraine susceptibility.")
md.append("- **FADS1**: Fatty acid desaturase genes adapted to agricultural diets (plant-rich vs meat-rich).")
md.append("- **ABCC11**: The variant responsible for dry earwax and reduced apocrine gland secretion, common in East Asia.")

md.append("## 3. High-Frequency and Differentiated Novel Candidates")
md.append(f"We scanned the entire dataset and identified **{len(candidates)}** candidate SNPs. Top 10 most differentiated:")
md.append(candidates.head(10).to_markdown(index=False))

md.append("## 4. Visualizations")
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

# By using \n\n, we ensure that markdown correctly closes tables before parsing the next element
md_text_for_md = "\n\n".join(md_for_md)
with open(MD_FILE, "w", encoding="utf-8") as f:
    f.write(md_text_for_md)

html_body = markdown.markdown("\n\n".join(md_for_html), extensions=['tables'])

html_full = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Adaptive Sweep Report</title>
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
        /* Scrollable tables */
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
        blockquote {{ border-left: 5px solid #3498db; padding: 10px 20px; color: #555; background: #f0f8ff; border-radius: 4px; }}
        img {{ max-width: 100%; height: auto; border-radius: 6px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin: 20px 0; display: block; }}
    </style>
</head>
<body>
<div class="container">
{html_body}
</div>
<script>
    // Wrap tables in a div to enable horizontal scrolling cleanly
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

print(f"Report regenerated with corrected spacing, gene descriptions, and JS-wrapped scrollable tables.")
