import json
import sys

pdb2domains = {}
for pdb in iter(input, ''):
    pdb2domains[pdb] = []

json.dump(pdb2domains, sys.stdout, indent=4)