import sys
import requests
import json
from os import path


list_file = '/home/adam/Workspace/C#/SecStrAnnot2_data/SecStrAPI/PF00067_rp15_accession_list.txt'
out_dir = '/home/adam/Workspace/C#/SecStrAnnot2_data/SecStrAPI/protein_api_cyps'
api_url = 'https://www.ebi.ac.uk/proteins/api/proteins/'

with open(list_file) as f:
    proteins = f.read().replace(',', '').split()

print(len(proteins), file = sys.stderr)

result = {}
failed = {}
for i, protein in enumerate(proteins):
    response = requests.get(api_url + protein)
    if response.status_code == 200:
        with open(path.join(out_dir, protein + '.json'), 'w') as w:
            w.write(response.text)
        print(i, protein, file=sys.stderr)
    else: 
        failed[protein] = response.status_code
        print(i, protein, 'failed', response.status_code, file=sys.stderr)

print('FAILED:', file=sys.stderr)
print(*failed.items(), sep='\n', file=sys.stderr)