'''
This Python3 script converts annotations in SecStrAPI format into tab-separated table.

Example usage:
    python3  annotation_json_to_tsv.py  annotations.json  --add_missing_sses
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
    row.extend(extra if extra is not None else 'NA' for extra in extras)
    if sse is not None:
        row.extend([sse[START], sse[END], length(sse), sse[TYPE]])
    else:
        row.extend(['NA', 'NA', 0, 'NA'])
    row.extend([longest_nested(sse, 'G'), longest_nested_starting(sse, 'G'), longest_nested_ending(sse, 'G'), 
        longest_nested(sse, 'H'), longest_nested(sse, 'I'), longest_nested_starting(sse, 'I'), longest_nested_ending(sse, 'I'), *count_GHI_bonds(sse)])
    return row

#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('all_annotations_file', help='JSON file with annotations in SecStrAPI format', type=str)
    parser.add_argument('--add_missing_sses', help='Add SSE labels which are missing for each domain, with all columns set to NA', action='store_true')
    parser.add_argument('--one_domain_per_pdb', help='Take only first domain for each PDB', action='store_true')
    args = parser.parse_args()
    return vars(args)


def main(all_annotations_file: str, add_missing_sses: bool = False, one_domain_per_pdb: bool = False) -> Optional[int]:
    '''Convert annotations in SecStrAPI format into tab-separated table.'''

    with open(all_annotations_file, 'r', encoding=lib.DEFAULT_ENCODING) as f:
        annotations = json.load(f)[ANNOTATIONS]

    labels = sorted(set(sse[LABEL] for pdb_annot in annotations.values() for dom_annot in pdb_annot.values() for sse in dom_annot[SSES]))

    result = []

    for pdb, pdb_annot in annotations.items():
        domain_annots = list(pdb_annot.items())
        if one_domain_per_pdb:
            domain_annots = domain_annots[:1]
        for domain, domain_annot in domain_annots:
            if add_missing_sses:
                for label in labels:
                    sse = next((sse for sse in domain_annot[SSES] if sse[LABEL] == label), None)
                    result.append(sse_to_row(sse, domain_annot[UNIPROT_ID], domain_annot[UNIPROT_NAME], pdb, domain, domain_annot[CHAIN], label))
            else:
                for sse in domain_annot[SSES]:
                    result.append(sse_to_row(sse, domain_annot[UNIPROT_ID], domain_annot[UNIPROT_NAME], pdb, domain, domain_annot[CHAIN], sse[LABEL]))
                    
    result.sort()

    print('UniProt', 'UniProt_name', 'PDB', 'Domain', CHAIN, LABEL, 'start', 'end', 'length', 'type', 
        'longest_G', 'longest_G_start', 'longest_G_end', 'longest_H', 'longest_I', 'longest_I_start', 'longest_I_end', 'bonds_G', 'bonds_H', 'bonds_I', 
        sep='\t')
    for row in result:
        print(*row, sep='\t')


if __name__ == '__main__':
    args = parse_args()
    exit_code = main(**args)
    if exit_code is not None:
        exit(exit_code)