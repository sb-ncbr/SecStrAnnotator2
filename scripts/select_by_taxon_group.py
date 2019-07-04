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
parser.add_argument('domain_list', help='JSON file with the list of domains in format (from merge_domain_lists.py)', type=str)
parser.add_argument('classification', help='TSV file with columns Domain, Taxon, Group', type=str)
parser.add_argument('group', help='Group to be selected', type=str)
args = parser.parse_args()

domain_list_file = args.domain_list
classification_file = args.classification
the_group = args.group

#  MAIN  ##############################################################################

with open(domain_list_file) as r:
    input_json = json.load(r)
pdb2domains = input_json[ANNOTATIONS]

domain2group = {}
with open(classification_file) as r:
    for line in r:
        line = line.strip()
        if not line.startswith('#'):
            domain, taxon, group = line.split()
            domain2group[domain] = group

for pdb, domains in list(pdb2domains.items()):
    for domain, annotation in list(domains.items()):
        if domain2group[domain] != the_group:
            domains.pop(domain)
    if len(domains) == 0:
        pdb2domains.pop(pdb)

print(json.dumps(input_json, indent=4))

n_pdbs = len(pdb2domains)
n_domains = sum( len(doms) for doms in pdb2domains.values() )
sys.stderr.write(f'Selected {n_domains} domains in {n_pdbs} PDB entries\n')