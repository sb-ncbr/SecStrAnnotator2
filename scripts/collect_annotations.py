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
INPUT_ALIGNMENT_EXT = '-alignment.json'
OUTPUT_EXT = '.json'
SEQUENCE_EXT = '.fasta'

USE_TWO_CLASS_SSE_TYPE = True
ADD_CONFIDENCE = True
METRIC_SIGNIFICANT_DIGITS = 2  # None for no rounding
DROP_FIELDS = ['start_vector', 'end_vector']

HIGH_CONFIDENCE_METRIC_THRESHOLD = 10.0
MEDIUM_CONFIDENCE_METRIC_THRESHOLD = 20.0

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('domain_list', help='JSON file with the list of domains in SecStrAPI format (from merge_domain_lists.py)', type=str)
parser.add_argument('input_directory', help='Directory with SSE annotations (from SecStrAnnotator2.dll)', type=str)
args = parser.parse_args()

domain_list_file = args.domain_list
input_directory = args.input_directory

#  FUNCTIONS  ###############################################################################

def two_class_sse_type(typ: str) -> str:
    if typ in 'hGHI':
        return 'h'
    elif typ in 'eEB':
        return 'e'
    else:
        raise Exception(f'Unknown SSE type {typ}')

def confidence_based_on_metric(metric_value: float) -> str:
    if metric_value <= HIGH_CONFIDENCE_METRIC_THRESHOLD:
        return 'high'
    elif metric_value <= MEDIUM_CONFIDENCE_METRIC_THRESHOLD:
        return 'medium'
    else:
        return 'low'

#  MAIN  ##############################################################################

with open(domain_list_file) as r:
    input_annotations = json.load(r)

api_version = input_annotations[API_VERSION]

for pdb, domains in input_annotations[ANNOTATIONS].items():
    convert_table_file = path.join(input_directory, pdb + '.label2auth.tsv')
    converter = lib.Label2AuthConverter(convert_table_file, unknown_ins_code_as_empty_string=True) if path.isfile(convert_table_file) else None
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
                if USE_TWO_CLASS_SSE_TYPE and 'type' in sse:
                    sse['type'] = two_class_sse_type(sse['type'])
                for field in DROP_FIELDS:
                    if field in sse:
                        sse.pop(field)
                if METRIC_SIGNIFICANT_DIGITS is not None and 'metric_value' in sse:
                    sse['metric_value'] = round(sse['metric_value'], METRIC_SIGNIFICANT_DIGITS)
                if ADD_CONFIDENCE and 'metric_value' in sse:
                    confidence = confidence_based_on_metric(sse['metric_value'])
                    lib.insert_after(sse, 'metric_value', [('confidence', confidence)])
        domain[CONNECTIVITY] = annot[CONNECTIVITY]
        try:
            with open(path.join(input_directory, name + INPUT_ALIGNMENT_EXT)) as r:
                domain[CANONICAL_ROTATION] = json.load(r)[CANONICAL_ROTATION]
        except IOError:
            pass
        if COMMENT in annot:
            domain[COMMENT] = annot[COMMENT]

json.dump(input_annotations, sys.stdout, indent=4)
print()

n_pdbs = len(input_annotations[ANNOTATIONS])
n_domains = sum( len(doms) for doms in input_annotations[ANNOTATIONS].values() )
sys.stderr.write(f'Collected annotations for {n_domains} domains in {n_pdbs} PDB entries\n')