import json
import csv
import os
import re

with open('/home/duri_bae/popgen_lab/raw/papers_fetched.json', 'r') as f:
    papers = json.load(f)

# Ensure at least 30 papers
if len(papers) < 30:
    print(f"Warning: Only found {len(papers)} papers. Expected at least 30.")

output_file = '/home/duri_bae/popgen_lab/output/AncientPaper.csv'
os.makedirs('/home/duri_bae/popgen_lab/output', exist_ok=True)

with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    # Header exactly matching the requirement plus DOI (implied by 9th column in example)
    writer.writerow(['Paper', 'Title', 'Summary', 'Problem Statement', 'Major Samples', 'Methodology', 'KeyFindings', 'Limitations&FutureWork', 'DOI'])
    
    count = 0
    for p in papers:
        if count >= 35:
            break
            
        authors = p.get('authors', '')
        # Basic cleanup of author name to avoid spaces/special chars in ID
        first_author = authors.split(',')[0].strip().split(' ')[-1] if authors else 'Unknown'
        # Remove non-alphanumeric from first author
        first_author = re.sub(r'[^a-zA-Z]', '', first_author)
        
        year = p.get('year', 'Unknown')
        journal = p.get('journal', 'Unknown')
        journal_clean = re.sub(r'[^a-zA-Z]', '', journal)
        paper_id = f"{first_author}{year}{journal_clean}"
        
        title = p.get('title', '')
        doi = p.get('doi', '')
        
        # Clean abstract HTML tags if any
        abstract = p.get('abstract', '')
        abstract = re.sub(r'<[^>]+>', '', abstract)
        
        if not abstract:
            summary = 'No abstract available.'
            prob = "unknown"
            findings = "needs_verification"
        else:
            # Splitting by sentences
            sentences = re.split(r'(?<=[.!?])\s+', abstract.strip())
            summary = "TL;DR: " + sentences[0]
            prob = sentences[1] if len(sentences) > 1 else "unknown"
            findings = " ".join(sentences[-2:]) if len(sentences) > 2 else abstract
            
        methodology = "Genome-wide ancient DNA sequencing analysis (inferred, needs_verification)."
        if "shotgun" in abstract.lower() or "capture" in abstract.lower():
            methodology = "Targeted capture or shotgun sequencing of ancient genomes."
            
        major_samples = "-"
        limitations = "needs_verification"
        
        writer.writerow([paper_id, title, summary, prob, major_samples, methodology, findings, limitations, doi])
        count += 1

print(f"Successfully processed {count} papers and saved to {output_file}")
