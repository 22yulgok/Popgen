import urllib.request
import urllib.parse
import json
import csv
import os

query = '"ancient DNA" AND ("Upper Paleolithic" OR "Neolithic" OR "hunter-gatherer" OR "Paleolithic" OR "LGM" OR "Last Glacial Maximum") AND (genome OR paleogenomics)'
encoded_query = urllib.parse.quote(query)
url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={encoded_query}&format=json&resultType=core&pageSize=40"

req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
try:
    with urllib.request.urlopen(req) as response:
        data = json.loads(response.read().decode('utf-8'))
except Exception as e:
    print(f"Error: {e}")
    exit(1)

papers = []
for result in data.get('resultList', {}).get('result', []):
    title = result.get('title', '')
    authors = result.get('authorString', '')
    year = result.get('pubYear', '')
    journal = result.get('journalTitle', '')
    doi = result.get('doi', '')
    abstract = result.get('abstractText', '')
    
    # Filter out empty DOIs
    if not doi:
        continue

    papers.append({
        'title': title,
        'authors': authors,
        'year': year,
        'journal': journal,
        'doi': doi,
        'abstract': abstract
    })

print(f"Fetched {len(papers)} papers.")

with open('/home/duri_bae/popgen_lab/raw/papers_fetched.json', 'w') as f:
    json.dump(papers[:35], f, indent=2)

print("Saved to /home/duri_bae/popgen_lab/raw/papers_fetched.json")
