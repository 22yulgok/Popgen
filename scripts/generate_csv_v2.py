import json
import csv
import os
import re

with open('/home/duri_bae/popgen_lab/raw/papers_fetched.json', 'r', encoding='utf-8') as f:
    papers = json.load(f)

output_file = '/home/duri_bae/popgen_lab/output/AncientPaper.csv'
os.makedirs('/home/duri_bae/popgen_lab/output', exist_ok=True)

# Use utf-8-sig to write the BOM so it perfectly matches the original file if it had one
with open(output_file, 'w', newline='', encoding='utf-8-sig') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Paper', 'Title', 'Summary', 'Problem Statement', 'Major Samples', 'Methodology', 'KeyFindings', 'Limitations&FutureWork', 'DOI'])
    
    count = 0
    for p in papers:
        if count >= 35:
            break
            
        authors = p.get('authors', '')
        # Extract first author's last name
        if authors:
            first_author_full = authors.split(',')[0].strip()
            first_author = first_author_full.split(' ')[-1]
            first_author = re.sub(r'[^a-zA-Z]', '', first_author)
        else:
            first_author = 'Unknown'
            
        year = p.get('year', 'Unknown')
        
        journal = p.get('journal', 'Unknown')
        # Title case and remove non-letters to make it like NatEcolEvol
        journal_clean = ''.join(word.capitalize() for word in re.sub(r'[^a-zA-Z\s]', '', journal).split())
        
        paper_id = f"{first_author}{year}{journal_clean}"
        
        title = p.get('title', '')
        doi = p.get('doi', '')
        
        abstract = p.get('abstract', '')
        abstract = re.sub(r'<[^>]+>', '', abstract).strip()
        
        if not abstract:
            summary = 'unknown'
            prob = "unknown"
            findings = "needs_verification"
        else:
            sentences = re.split(r'(?<=[.!?])\s+', abstract)
            summary = sentences[0] if sentences else "unknown"
            prob = sentences[1] if len(sentences) > 1 else "unknown"
            findings = " ".join(sentences[-2:]) if len(sentences) > 2 else "needs_verification"
            
        methodology = "needs_verification"
        if "shotgun" in abstract.lower() or "capture" in abstract.lower() or "genome" in abstract.lower() or "dna" in abstract.lower():
            methodology = "Genome-wide ancient DNA sequencing analysis (inferred, needs_verification)"
            
        major_samples = "-"
        limitations = "needs_verification"
        
        writer.writerow([paper_id, title, summary, prob, major_samples, methodology, findings, limitations, doi])
        count += 1

print(f"Successfully processed {count} papers and saved to {output_file}")
