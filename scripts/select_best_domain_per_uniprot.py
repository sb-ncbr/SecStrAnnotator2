import os
import sys
import json
import requests
import argparse
from collections import defaultdict
from typing import Tuple
from os import path

import lib

#  CONSTANTS  ##############################################################################

from constants import *

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('domain_list', help='JSON file with the list of domains in SecStrAPI format (from merge_domain_lists.py)', type=str)
args = parser.parse_args()

domain_list_file = args.domain_list

#  MAIN  ##############################################################################

with open(domain_list_file) as r:
    input_json = json.load(r)

pdb2domains = input_json[ANNOTATIONS]

pdbs = pdb2domains.keys()
pdb2quality = { pdb: lib.get_structure_quality(pdb) for pdb in pdbs }

DOMAINS_IN_DICT = len(pdb2domains) > 0 and isinstance(next(iter(pdb2domains.values())), dict)

uniprot2domains = defaultdict(list)
for pdb, domains in pdb2domains.items():
    if DOMAINS_IN_DICT:
        pdbs_domains_keys = [ (pdb, domain, name) for name, domain in domains.items() ]
    else:
        pdbs_domains_keys = [ (pdb, domain) for domain in domains ]
    for pdb_dom_key in pdbs_domains_keys:
        uniprot = pdb_dom_key[1][UNIPROT_ID]
        uniprot2domains[uniprot].append(pdb_dom_key)

pdb2best_domains = defaultdict(dict) if DOMAINS_IN_DICT else defaultdict(list)

for pdbs_domains_keys in uniprot2domains.values():
    best = max(pdbs_domains_keys, key=lambda pdb_dom_key: pdb2quality[pdb_dom_key[0]])
    if DOMAINS_IN_DICT:
        pdb, domain, key = best
        pdb2best_domains[pdb][key] = domain
    else:
        pdb, domain = best
        pdb2best_domains[pdb].append(domain)

input_json[ANNOTATIONS] = { pdb: domains for (pdb, domains) in sorted(pdb2best_domains.items()) }

json.dump(input_json, sys.stdout, indent=4)
print()

n_pdbs = len(pdb2best_domains)
n_domains = sum( len(doms) for doms in pdb2best_domains.values() )
sys.stderr.write(f'Selected {n_domains} domains in {n_pdbs} PDB entries\n')