'''
This Python3 script removes unwanted fields from annotation file in SecStrAPI format.

Example usage:
    python3  annotations.json
'''

import argparse
from typing import Dict, Any, Optional, Union, List
import os
import sys
import json
from os import path
import shutil

import lib
from constants import *

#  CONSTANTS  ################################################################################

DEFAULT_FIELDS = 'label,chain_id,start,end,auth_chain_id,auth_start,auth_start_ins_code,auth_end,auth_end_ins_code,type,color,sheet_id,metric_value,confidence,sequence,reference_residue'
OUTPUT_EXT = '.json'

#  FUNCTIONS  ################################################################################


#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_annotations_file', help='JSON file with the list of domains with annotations in SecStrAPI format', type=str)
    parser.add_argument('--fields', help=f'Comma-separated list of fields which should be kept in "secondary_structure_elements" (default: {DEFAULT_FIELDS})', type=str, default=DEFAULT_FIELDS)
    args = parser.parse_args()
    return vars(args)


def main(input_annotations_file: str, fields: Union[str, List[str]]=DEFAULT_FIELDS) -> Optional[int]:
    '''Remove unwanted fields from annotation file in SecStrAPI format (keep only fields listed in parameter fields).'''
    if isinstance(fields, str):
        fields = fields.split(',')
    
    with open(input_annotations_file, 'r', encoding=lib.DEFAULT_ENCODING) as r:
        all_annotations = json.load(r)

    pdb2domains = all_annotations[ANNOTATIONS]
    for pdb, domains in pdb2domains.items():
        for domain, annot in domains.items():
            for sse in annot[SSES]:
                for field in list(sse.keys()):
                    if field not in fields:
                        sse.pop(field)

    json.dump(all_annotations, sys.stdout, indent=4)
    print()


if __name__ == '__main__':
    args = parse_args()
    exit_code = main(**args)
    if exit_code is not None:
        exit(exit_code)