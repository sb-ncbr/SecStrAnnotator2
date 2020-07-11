'''
This Python3 script reads domains in SecStrAPI format, selects only domains belonging to specified taxonomy group, and prints them to output in SecStrAPI format.

Example usage:
    python3  select_by_taxon_group.py  domains.json  domains_taxons_groups.tsv  Bact
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
    parser.add_argument('domain_list_file', help='JSON file with the list of domains in format (from merge_domain_lists.py)', type=str)
    parser.add_argument('classification_file', help='TSV file with columns Domain, Taxon, Group', type=str)
    parser.add_argument('the_group', help='Group to be selected', type=str)
    args = parser.parse_args()
    return vars(args)


def main(domain_list_file: str, classification_file: str, the_group: str) -> None:
    '''Read domains in SecStrAPI format, select only domains belonging to specified taxonomy group, and print them to output in SecStrAPI format.'''
    with open(domain_list_file, 'r', encoding=lib.DEFAULT_ENCODING) as r:
        input_json = json.load(r)
    pdb2domains = input_json[ANNOTATIONS]

    domain2group = {}
    with open(classification_file, 'r', encoding=lib.DEFAULT_ENCODING) as r:
        for line in r:
            line = line.strip()
            if not line.startswith('#'):
                domain, taxon, group = line.split()
                domain2group[domain] = group

    for pdb, domains in list(pdb2domains.items()):
        for domain, annotation in list(domains.items()):
            if domain2group[domain] != the_group:
                domains.pop(domain)
        if len(domains) == 0:
            pdb2domains.pop(pdb)

    print(json.dumps(input_json, indent=4))

    n_pdbs = len(pdb2domains)
    n_domains = sum( len(doms) for doms in pdb2domains.values() )
    sys.stderr.write(f'Selected {n_domains} domains in {n_pdbs} PDB entries\n')


if __name__ == '__main__':
    args = parse_args()
    main(**args)