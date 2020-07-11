'''
This Python3 script reads domains in SecStrAPI format and prints them in simple JSON format { pdb: [[domain_name, chain, range]] }

Example usage:
    python3  simplify_domain_list.py  domains.json
'''

import argparse
from typing import Dict, Any
import sys
import json

import lib
from constants import *

#  CONSTANTS  ################################################################################


#  FUNCTIONS  ################################################################################


#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_file', help='Domain list in SecStrAPI format', type=str)
    args = parser.parse_args()
    return vars(args)


def main(input_file: str) -> None:
    '''Read domains in SecStrAPI format from input_file and print them in simple JSON format { pdb: [[domain_name, chain, range]] }'''
    with open(input_file, 'r', encoding=lib.DEFAULT_ENCODING) as r:
        domain_list = json.load(r)
        pdb2domains = domain_list[ANNOTATIONS]

    simple_list = {}
    for pdb, domains in pdb2domains.items():
        simple_list[pdb] = [ (name, dom[CHAIN], dom[RANGES]) for name, dom in lib.iterate_names_domains(domains) ] 

    json.dump(simple_list, sys.stdout, indent=4)
    print()

    n_pdbs = len(simple_list)
    n_domains = sum( len(doms) for doms in simple_list.values() )
    sys.stderr.write(f'Formatted {n_domains} domains in {n_pdbs} PDB entries\n')


if __name__ == '__main__':
    args = parse_args()
    main(**args)