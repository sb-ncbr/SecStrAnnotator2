'''
This Python3 script reads domains in SecStrAPI format, selects one best-quality domain for each UniProtID, and prints them to output in SecStrAPI format.

Example usage:
    python3  select_best_domain_per_uniprot.py  domains.json 
'''

import argparse
from typing import Dict, Any, Tuple
import os
from os import path
import sys
import json
import requests
from collections import defaultdict

import lib
from constants import *

#  CONSTANTS  ################################################################################


#  FUNCTIONS  ################################################################################


#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('domain_list_file', help='JSON file with the list of domains in SecStrAPI format (from domain_lists_to_SecStrAPI_format.py)', type=str)
    args = parser.parse_args()
    return vars(args)


def main(domain_list_file: str) -> None:
    '''Read domains in SecStrAPI format, select one best-quality domain for each UniProtID, and print them to output in SecStrAPI format.'''
    with open(domain_list_file, 'r', encoding=lib.DEFAULT_ENCODING) as r:
        input_json = json.load(r)

    pdb2domains = input_json[ANNOTATIONS]

    pdbs = pdb2domains.keys()
    pdb2quality = {}
    progress_bar = lib.ProgressBar(len(pdbs), title=f'Getting structure quality for {len(pdbs)} PDB entries', writer=sys.stderr).start()
    for pdb in pdbs:
        pdb2quality[pdb] = lib.get_structure_quality(pdb)
        progress_bar.step()
    progress_bar.finalize()

    DOMAINS_IN_DICT = len(pdb2domains) > 0 and isinstance(next(iter(pdb2domains.values())), dict)

    uniprot2domains = defaultdict(list)
    for pdb, domains in pdb2domains.items():
        if DOMAINS_IN_DICT:
            pdbs_domains_keys = [ (pdb, domain, name) for name, domain in domains.items() ]
        else:
            pdbs_domains_keys = [ (pdb, domain) for domain in domains ]
        for pdb_dom_key in pdbs_domains_keys:
            uniprot = pdb_dom_key[1][UNIPROT_ID]
            if uniprot is not None:
                uniprot2domains[uniprot].append(pdb_dom_key)

    pdb2best_domains = defaultdict(dict) if DOMAINS_IN_DICT else defaultdict(list)

    for pdbs_domains_keys in uniprot2domains.values():
        best = max(pdbs_domains_keys, key=lambda pdb_dom_key: pdb2quality[pdb_dom_key[0]])
        if DOMAINS_IN_DICT:
            pdb, domain, key = best
            pdb2best_domains[pdb][key] = domain
        else:
            pdb, domain = best
            pdb2best_domains[pdb].append(domain)
        # sys.stderr.write(f'{pdb}\t{pdb2quality[pdb]}\n')

    input_json[ANNOTATIONS] = { pdb: domains for (pdb, domains) in sorted(pdb2best_domains.items()) }

    json.dump(input_json, sys.stdout, indent=4)
    print()

    n_pdbs = len(pdb2best_domains)
    n_domains = sum( len(doms) for doms in pdb2best_domains.values() )
    sys.stderr.write(f'Selected {n_domains} domains in {n_pdbs} PDB entries\n')


if __name__ == '__main__':
    args = parse_args()
    main(**args)