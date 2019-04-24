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

INPUT_EXT = '-annotated.sses.json'
OUTPUT_EXT = '.json'
SEQUENCE_EXT = '.fasta'

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

for pdb, domains in input_annotations[ANNOTATIONS].items():
    convert_table_file = path.join(input_directory, pdb + '.label2auth.tsv')
    converter = lib.Label2AuthConverter(convert_table_file, unknown_ins_code_as_empty_string=True) if path.isfile(convert_table_file) else None
    # dom_list = domains.values() if isinstance(domains, dict) else domains[:]
    # name_domain_gen = domains.items() if isinstance(domains, dict) else ( (f'{dom[PDB]},{dom[CHAIN]},{dom[RANGES]}', dom) for dom in domains )
    for name, domain in lib.iterate_names_domains(domains):
        with open(path.join(input_directory, name + INPUT_EXT)) as r:
            annot = json.load(r)[pdb]
        if converter is not None:
            auth_chain, auth_ranges = converter.auth_chain_ranges(domain[CHAIN], domain[RANGES])
            lib.insert_after(domain, RANGES, ((AUTH_CHAIN, auth_chain), (AUTH_RANGES, auth_ranges)))
        domain[SSES] = annot[SSES]
        if converter is not None:
            for sse in domain[SSES]:
                auth_chain, auth_start, auth_start_ins = converter.auth_chain_resi_ins(sse[CHAIN_ID], sse[START])
                auth_chain, auth_end, auth_end_ins = converter.auth_chain_resi_ins(sse[CHAIN_ID], sse[END])
                lib.insert_after(sse, END, 
                    ((AUTH_CHAIN_ID, auth_chain), (AUTH_START, auth_start), (AUTH_START_INS, auth_start_ins), (AUTH_END, auth_end), (AUTH_END_INS, auth_end_ins)))
        domain[CONNECTIVITY] = annot[CONNECTIVITY]
        if COMMENT in annot:
            domain[COMMENT] = annot[COMMENT]

json.dump(input_annotations, sys.stdout, indent=4)
print()

n_pdbs = len(input_annotations[ANNOTATIONS])
n_domains = sum( len(doms) for doms in input_annotations[ANNOTATIONS].values() )
sys.stderr.write(f'Collected annotations for {n_domains} domains in {n_pdbs} PDB entries\n')