'''
This Python3 script takes SSE annotations from a single file in SecStrAPI format and creates separate file for each PDB entry.

Example usage:
    python3  divide_annotations_by_pdb.py  annotations.json  annotations/  --min_dir annotations_min/
'''

import argparse
from typing import Dict, Any, Optional
import os
import sys
import json
from os import path
import shutil

import lib
from constants import *

#  CONSTANTS  ################################################################################

OUTPUT_EXT = '.json'

#  FUNCTIONS  ################################################################################


#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('all_annotations_file', help='JSON file with the list of domains with annotations in SecStrAPI format', type=str)
    parser.add_argument('output_directory', help='Directory to output per-PDB-entry annotations', type=str)
    parser.add_argument('--min_dir', help='Directory to output per-PDB-entry annotations in minified JSON', type=str, default=None)
    args = parser.parse_args()
    return vars(args)


def main(all_annotations_file: str, output_directory: str, min_dir: Optional[str] = None) -> Optional[int]:
    '''Takes SSE annotations from a single file in SecStrAPI format and creates separate file for each PDB entry.'''

    with open(all_annotations_file, 'r', encoding=lib.DEFAULT_ENCODING) as r:
        all_annotations = json.load(r)

    shutil.rmtree(output_directory, ignore_errors=True)
    os.makedirs(output_directory)
    if min_dir is not None:
        shutil.rmtree(min_dir, ignore_errors=True)
        os.makedirs(min_dir)

    api_version = all_annotations[API_VERSION]
    pdb2domains = all_annotations[ANNOTATIONS]
    for pdb, domains in pdb2domains.items():
        result = { API_VERSION: api_version, ANNOTATIONS: { pdb: domains } }
        with open(path.join(output_directory, pdb + OUTPUT_EXT), 'w', encoding=lib.DEFAULT_ENCODING) as w:
            json.dump(result, w, indent=4)
            w.write('\n')
        if min_dir is not None:
            with open(path.join(min_dir, pdb + OUTPUT_EXT), 'w', encoding=lib.DEFAULT_ENCODING) as w:
                json.dump(result, w, separators=(',', ':'))
                w.write('\n')

    n_pdbs = len(all_annotations[ANNOTATIONS])
    n_domains = sum( len(doms) for doms in all_annotations[ANNOTATIONS].values() )
    sys.stderr.write(f'Divided annotations for {n_domains} domains in {n_pdbs} PDB entries\n')


if __name__ == '__main__':
    args = parse_args()
    exit_code = main(**args)
    if exit_code is not None:
        exit(exit_code)