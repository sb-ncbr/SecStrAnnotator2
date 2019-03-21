import os
import sys
import json
import requests
import argparse
from collections import defaultdict
from typing import Tuple
from os import path

#  CONSTANTS  ##############################################################################

from constants import *

INPUT_EXT = '-annotated.sses.json'
OUTPUT_EXT = '.json'
SEQUENCE_EXT = '.fasta'

#  FUNCTIONS  ##############################################################################

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('domain_list', help='JSON file with the list of domains in SecStrAPI format (from merge_domain_lists.py)', type=str)
parser.add_argument('input_directory', help='Directory with SSE annotations (from SecStrAnnot2.dll)', type=str)
args = parser.parse_args()

domain_list_file = args.domain_list
input_directory = args.input_directory

#  MAIN  ##############################################################################

with open(domain_list_file) as r:
    input_annotations = json.load(r)

api_version = input_annotations[API_VERSION]

label2domain_sequence = defaultdict(lambda: [])

for pdb, domains in input_annotations[ANNOTATIONS].items():
    dom_list = domains.values() if isinstance(domains, dict) else domains[:]
    for domain in dom_list:
        name = ','.join((domain[PDB], domain[CHAIN], domain[RANGES]))
        with open(path.join(input_directory, name + INPUT_EXT)) as r:
            annot = json.load(r)[pdb]
        domain[SSES] = annot[SSES]
        domain[CONNECTIVITY] = annot[CONNECTIVITY]
        if COMMENT in annot:
            domain[COMMENT] = annot[COMMENT]
        for sse in domain[SSES]:
            if LABEL in sse and sse[LABEL] is not None:
                label2domain_sequence[sse[LABEL]].append((name, sse[SEQUENCE]))

json.dump(input_annotations, sys.stdout, indent=4)
print()

n_pdbs = len(input_annotations[ANNOTATIONS])
n_domains = sum( len(doms) for doms in input_annotations[ANNOTATIONS].values() )
sys.stderr.write(f'Collected annotations for {n_domains} domains in {n_pdbs} PDB entries\n')