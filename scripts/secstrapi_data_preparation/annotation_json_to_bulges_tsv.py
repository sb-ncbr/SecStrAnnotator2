'''
This Python3 script extracts beta-bulges from annotations in SecStrAPI format and prints them in tab-separated table.

Example usage:
    python3  annotation_json_to_tsv.py  annotations.json
'''

import argparse
from typing import Dict, Any, Optional
import json
import os
import sys
import math

import lib
from constants import *

#  CONSTANTS  ################################################################################

BULGE_TYPES='nNmMtTsSoOpPqQrRlL'

#  FUNCTIONS  ################################################################################

def equal(sse1, sse2):
    return difference(sse1, sse2) == 0
    if sse1 == None or sse2 == None:
        return False
    else:
        return sse1[START] == sse2[START] and sse1[END] == sse2[END]

def difference(sse1, sse2):
    if sse1 == None or sse2 == None:
        return math.nan
    else:
        return abs(sse1[START] - sse2[START]) + abs(sse1[END] - sse2[END])

def length(sse):
    return sse[END] - sse[START] + 1

def longest_nested(sse, nested_type):
    if sse is not None and NESTED_SSES in sse:
        nested = sse[NESTED_SSES]
        lengths = [length(nes) for nes in nested if nes[TYPE]==nested_type]
        # sys.stderr.write('Lengths: ' + str(lengths) + '\n')
        longest = max(lengths + [0])
        return longest
    else:
        return 0

def longest_nested_starting(sse, nested_type, from_start_tolerance=1):
    if sse is not None and NESTED_SSES in sse:
        nested = sse[NESTED_SSES]
        lengths = [length(nes) for nes in nested if nes[TYPE]==nested_type and nes[START]<=sse[START]+from_start_tolerance]
        # sys.stderr.write('Lengths: ' + str(lengths) + '\n')
        longest = max(lengths + [0])
        return longest
    else:
        return 0

def longest_nested_ending(sse, nested_type, from_end_tolerance=1):
    if sse is not None and NESTED_SSES in sse:
        nested = sse[NESTED_SSES]
        lengths = [length(nes) for nes in nested if nes[TYPE]==nested_type and nes[END]>=sse[END]-from_end_tolerance]
        # sys.stderr.write('Lengths: ' + str(lengths) + '\n')
        longest = max(lengths + [0])
        return longest
    else:
        return 0

def elemwise_sum(lists):
    return [ sum(elems) for elems in zip(*lists) ]

def count_GHI_bonds(sse):
    if sse is None:
        return [0, 0, 0]
    elif NESTED_SSES in sse:
        return elemwise_sum( count_GHI_bonds(nested) for nested in sse[NESTED_SSES] )
    else:
        if sse[TYPE] == 'G':
            return [length(sse) - 1, 0, 0]
        elif sse[TYPE] == 'H':
            return [0, length(sse) - 2, 0]
        elif sse[TYPE] == 'I':
            return [0, 0, length(sse) - 3]
        else:
            return [0, 0, 0]

def sse_to_row(sse, *extras):
    row = []
    row.extend( extra if extra is not None else 'NA' for extra in extras )
    row.extend( sse[field] if sse is not None else 'NA' for field in fields )
    row.extend([longest_nested(sse, 'G'), longest_nested_starting(sse, 'G'), longest_nested_ending(sse, 'G'), 
        longest_nested(sse, 'H'), longest_nested(sse, 'I'), longest_nested_starting(sse, 'I'), longest_nested_ending(sse, 'I'), *count_GHI_bonds(sse)])
    return row

def find_bulges_recursive(sse):
    result = []
    if sse[TYPE] in BULGE_TYPES:
        result.append(sse)
    for nested in sse.get(NESTED_SSES, []):
        result.extend(find_bulges_recursive(nested))
    return result

#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('all_annotations_file', help='JSON file with annotations in SecStrAPI format', type=str)
    parser.add_argument('--one_domain_per_pdb', help='Take only first domain for each PDB', action='store_true')
    args = parser.parse_args()
    return vars(args)


def main(all_annotations_file: str, one_domain_per_pdb: bool = False) -> Optional[int]:
    '''Extract beta-bulges from annotations in SecStrAPI format and print them in tab-separated table.'''
    
    with open(all_annotations_file, 'r', encoding=lib.DEFAULT_ENCODING) as f:
        annotations=json.load(f)[ANNOTATIONS]

    labels = sorted(set( sse[LABEL] for pdb_annot in annotations.values() for dom_annot in pdb_annot.values() for sse in dom_annot[SSES] ))
    result = []

    for pdb, pdb_annot in annotations.items():
        domain_annots = list(pdb_annot.items())
        if one_domain_per_pdb:
            domain_annots = domain_annots[:1]
        for domain, domain_annot in domain_annots:
            uni = domain_annot[UNIPROT_ID]
            sse_bulges = [ (sse, bulge) for sse in domain_annot[SSES] for bulge in find_bulges_recursive(sse) ]
            for sse, bulge in sse_bulges:
                row = [uni, pdb, domain, sse[LABEL], bulge[LABEL], bulge[CHAIN], bulge[START], bulge[END], length(bulge), bulge[TYPE]]
                result.append(row)
                    
    print('UniProt', 'PDB', 'Domain', 'label', 'bulge_label', CHAIN, 'start', 'end', 'length', 'type', sep='\t')
    for row in sorted(result):
        print(*row, sep='\t')


if __name__ == '__main__':
    args = parse_args()
    exit_code = main(**args)
    if exit_code is not None:
        exit(exit_code)