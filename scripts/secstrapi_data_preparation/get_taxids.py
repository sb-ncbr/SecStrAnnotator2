'''
This Python3 script find taxonomy ID for each protein domain in the input file and prints a tab-separated table domain|taxid.
'''

import argparse
from typing import Dict, Any
import os
import sys
import json
import requests

import lib
from constants import *

#  CONSTANTS  ################################################################################

URL = 'http://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/'

#  FUNCTIONS  ################################################################################


#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('domain_list_file', help='JSON file with the list of domains in SecStrAPI format (from merge_domain_lists.py)', type=str)
    args = parser.parse_args()
    return vars(args)


def main(domain_list_file: str) -> None:
    '''Find taxonomy ID for each protein domain in domain_list_file and print a tab-separated table domain|taxid.'''
    with open(domain_list_file, 'r', encoding=lib.DEFAULT_ENCODING) as r:
        input_annotations = json.load(r)[ANNOTATIONS]

    progress_bar = lib.ProgressBar(len(input_annotations), title=f'Getting taxonomy ID for {len(input_annotations)} PDB entries', writer=sys.stderr).start()
    for pdb, pdb_annot in input_annotations.items():
        response = json.loads(requests.get(URL + pdb).text)
        for domain, domain_annot in pdb_annot.items():
            chain = domain_annot[CHAIN]
            entities = [ entity for entity in response[pdb] if chain in entity['in_struct_asyms'] ] 
            taxids = [ source['tax_id'] for entity in entities for source in entity.get('source',[]) ]
            if len(set(taxids))!=1:
                sys.stderr.write('WARNING: ' + pdb + ' has non-unique taxon ID: ' + ', '.join([str(t) for t in taxids]) + '\n')
            print(domain, taxids[0], sep='\t')
        progress_bar.step()
    progress_bar.finalize()


if __name__ == '__main__':
    args = parse_args()
    main(**args)