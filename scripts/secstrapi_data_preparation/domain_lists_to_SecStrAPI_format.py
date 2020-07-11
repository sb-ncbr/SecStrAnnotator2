'''
This Python3 script merges domains from multiple files in simple JSON format { pdb: [[domain_name, chain, range]] }
and prints them in SecStrAPI format.

Example usage:
    python3  domain_lists_to_SecStrAPI_format  CATH  1.10.630.10  set_cath.simple.json  Pfam  PF00067  ste_pfam.simple.json  --api_version 1.0
'''

import argparse
from typing import Dict, Any, Tuple, List
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
    parser.add_argument('lists', help='One or more triples of arguments SOURCE FAMILY FILE, where SOURCE is the name of a source (CATH, Pfam...), FAMILY is the identifier of the family in the source (1.10.630.10, PF00067...), and FILE is a JSON file with the list of domains from the source in format {PDB:[[domain, chain, ranges]]}', nargs='*')
    parser.add_argument('--api_version', help='API version information to include in the output', type=str, default=None)
    parser.add_argument('--domain_naming', help='Type of chain ID to put into the created domain name', choices=['label', 'auth'], default='auth')
    args = parser.parse_args()
    return vars(args)


def main(lists: List[str], api_version: str = None, domain_naming: str = 'auth') -> None:
    '''Merge domains from multiple files in simple JSON format { pdb: [[domain_name, chain, range]] } and print them in SecStrAPI format'''

    if len(lists) % 3 != 0:
        parser.error('The number of LISTS arguments must be a multiple of 3')
    else:
        sources = lists[0::3]
        families = lists[1::3]
        files = lists[2::3]
    if api_version is None:
        sys.stderr.write('WARNING: "api_version" is set to null\n')

    annotations = {}
    uniprot_manager = lib.UniProtManager()

    pdb2domains_dicts = []
    pdbs_todo = set()
    for filename in files:
        with open(filename, 'r', encoding=lib.DEFAULT_ENCODING) as r:
            pdb2domains = json.load(r)
            pdb2domains_dicts.append(pdb2domains)
            pdbs_todo.update(pdb2domains.keys())

    progress_bar = lib.ProgressBar(len(pdbs_todo), title=f'Formatting {len(pdbs_todo)} PDBs and getting UniProt info', writer=sys.stderr).start()
    for source, family, pdb2domains in zip(sources, families, pdb2domains_dicts):
        for pdb, domains in pdb2domains.items():
            if pdb not in annotations:
                annotations[pdb] = {} if DOMAINS_IN_DICT else []
            for domain in domains:
                domain_name = domain[DOMAIN_NAME]
                chain = domain[CHAIN]
                ranges = domain[RANGES]
                auth_chain = domain[AUTH_CHAIN]
                auth_ranges = domain[AUTH_RANGES]
                if DOMAINS_IN_DICT:
                    to_unify = [ dom for dom in annotations[pdb].values() if dom[CHAIN] == chain ]
                else:
                    to_unify = [ dom for dom in annotations[pdb] if dom[CHAIN] == chain ]
                if len(to_unify) == 0:  # create a new domain
                    our_domain_id = lib.create_domain_id(pdb, chain if domain_naming=='label' else auth_chain)
                    dom = { PDB: pdb, CHAIN: chain, RANGES: ':', AUTH_CHAIN: auth_chain, AUTH_RANGES: ':' }  # Take whole chain always
                    dom[UNIPROT_ID], dom[UNIPROT_NAME] = uniprot_manager.get_uniprot_id_and_name(pdb, chain)
                    dom[MAPPINGS] = []
                    if DOMAINS_IN_DICT:
                        annotations[pdb][our_domain_id] = dom
                    else:
                        annotations[pdb].append(dom)
                elif len(to_unify) == 1:  # just add a mapping to an existing domain
                    dom = to_unify[0]
                else:
                    raise Exception(f'Domain {domain_name} ({pdb} {chain} {ranges}) can be unified with more than one other domain')
                new_mapping = {DOMAIN_NAME: domain_name, SOURCE: source, FAMILY_ID: family, PDB: pdb, CHAIN: chain, RANGES: ranges, AUTH_CHAIN: auth_chain, AUTH_RANGES: auth_ranges}
                dom[MAPPINGS].append(new_mapping)
            if pdb in pdbs_todo:
                progress_bar.step()
                pdbs_todo.discard(pdb)
    progress_bar.finalize()

    annotations = { pdb: annot for pdb, annot in sorted(annotations.items()) }

    result = { API_VERSION: api_version, ANNOTATIONS: annotations }
    json.dump(result, sys.stdout, indent=4)
    print()

    n_pdbs = len(annotations)
    n_domains = sum( len(doms) for doms in annotations.values() )
    sys.stderr.write(f'Formatted {n_domains} domains in {n_pdbs} PDB entries\n')


if __name__ == '__main__':
    args = parse_args()
    main(**args)