# Example of usage:
# python3 ../../SecStrAnnot2/scripts/merge_domain_lists.py CATH 1.10.630.10 cyps_cath_20190312.json Pfam PF00067 cyps_pfam_20190312.json --api_version 1.0 > cyps_merged_20190312.json

import os
import sys
import json
import requests
import argparse
from collections import defaultdict
from typing import Tuple
from os import path

#  CONSTANTS  ##############################################################################

INPUT_EXT = '-annotated.sses.json'
OUTPUT_EXT = '.json'
SEQUENCE_EXT = '.fasta'

API_VERSION = 'api_version'
ANNOTATIONS = 'annotations'
PDB = 'pdb'
CHAIN = 'chain'
RANGES = 'ranges'
UNIPROT_ID = 'uniprot_id'
UNIPROT_NAME = 'uniprot_name'
MAPPINGS = 'domain_mappings'
DOMAIN_NAME = 'domain'
SOURCE = 'source'
FAMILY_ID = 'family'

SSES = 'secondary_structure_elements'
CONNECTIVITY = 'beta_connectivity'
COMMENT = 'comment'
LABEL = 'label'
SEQUENCE = 'sequence'

#  FUNCTIONS  ##############################################################################

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('domain_list', help='JSON file with the list of domains in SecStrAPI format (from merge_domain_lists.py)', type=str)
parser.add_argument('input_directory', help='Directory with SSE annotations (from SecStrAnnot2.dll)', type=str)
parser.add_argument('output_directory', help='Directory for SSE annotations in SecStrAPI format', type=str)
parser.add_argument('--sequences_directory', help='Directory for SSE sequences in multi-FASTA format', type=str, default=None)
args = parser.parse_args()

domain_list_file = args.domain_list
input_directory = args.input_directory
output_directory = args.output_directory
sequences_directory = args.sequences_directory

#  MAIN  ##############################################################################

with open(domain_list_file) as r:
    input_annotations = json.load(r)

os.makedirs(output_directory, exist_ok=True)
if sequences_directory is not None:
    os.makedirs(sequences_directory, exist_ok=True)

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
    with open(path.join(output_directory, pdb + OUTPUT_EXT), 'w') as w:
        result = { API_VERSION: api_version, ANNOTATIONS: { pdb: domains } }
        json.dump(result, w, indent=4)

with open(path.join(output_directory, 'all' + OUTPUT_EXT), 'w') as w:
    json.dump(input_annotations, w, indent=4)

if sequences_directory is not None:
    for label, domnames_sequences in label2domain_sequence.items():
        with open(path.join(sequences_directory, label + SEQUENCE_EXT), 'w') as w:
            for domname, sequence in domnames_sequences:
                w.write(f'>{domname}\n{sequence}\n')
