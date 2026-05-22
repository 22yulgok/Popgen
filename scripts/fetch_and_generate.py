import urllib.request
import urllib.parse
import json
import csv
import os
import re
import time

query = '"ancient DNA" AND ("Upper Paleolithic" OR "Neolithic" OR "hunter-gatherer" OR "Paleolithic" OR "LGM" OR "Last Glacial Maximum") AND (genome OR paleogenomics)'
encoded_query = urllib.parse.quote(query)
url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={encoded_query}&format=json&resultType=core&pageSize=40"

req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    with urllib.request.urlopen(req) as response:
        data = json.loads(response.read().decode('utf-8'))
except Exception as e:
    print(f"Error fetching data: {e}")
    exit(1)

papers = data.get('resultList', {}).get('result', [])
if len(papers) < 30:
    print(f"Warning: Only found {len(papers)} papers. Expected at least 30.")

output_file = '/home/duri_bae/popgen_lab/output/AncientPaper.csv'
os.makedirs('/home/duri_bae/popgen_lab/output', exist_ok=True)

with open(output_file, 'w', newline='', encoding='utf-8-sig') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Paper', 'Title', 'Summary', 'Problem Statement', 'Major Samples', 'Methodology', 'KeyFindings', 'Limitations&FutureWork', 'DOI'])
    
    count = 0
    for p in papers:
        if count >= 35:
            break
            
        authors = p.get('authorString', '')
        if authors:
            first_author_full = authors.split(',')[0].strip()
            # In EuropePMC, authorString is usually "Lastname Initial, ...". The last name is the first word.
            first_author = first_author_full.split(' ')[0]
            first_author = re.sub(r'[^a-zA-Z]', '', first_author)
        else:
            first_author = 'Unknown'
            
        year = str(p.get('pubYear', 'Unknown'))
        
        journal = p.get('journalTitle', '')
        if not journal:
            # try finding it in journalInfo
            journalInfo = p.get('journalInfo', {})
            if isinstance(journalInfo, dict):
                journal = journalInfo.get('journal', {}).get('title', '')
        if not journal:
            journal = 'Unknown'
            
        if journal == 'Unknown':
            journal_clean = 'Unknown'
        else:
            journal_clean = ''.join(word.capitalize() for word in re.sub(r'[^a-zA-Z\s]', '', journal).split())
        
        paper_id = f"{first_author}{year}{journal_clean}"
        
        title = p.get('title', '')
        doi = p.get('doi', '')
        if not doi:
            continue
            
        abstract = p.get('abstractText', '')
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
