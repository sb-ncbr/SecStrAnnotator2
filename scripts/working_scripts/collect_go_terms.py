'''
This Python3 script does foo ...

Example usage:
    python3  foo.py  --foo 4  foo.txt 
'''
# TODO add description and example usage in docstring

import sys
import argparse
import json
import requests
from collections import defaultdict
from typing import Dict, Any, Optional

#  CONSTANTS  ################################################################################

API_URL = 'https://www.ebi.ac.uk/pdbe/api/mappings/{pdb}'

#  FUNCTIONS  ################################################################################


#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_json', help='JSON file with input PDBs and UniProts (like set_NR.json)', type=str)
    args = parser.parse_args()
    return vars(args)


def main(input_json: str) -> Optional[int]:
    # TODO add parameters
    '''Foo'''
    # TODO add docstring

    # Read domains from input file
    with open(input_json) as r:
        inp = json.load(r)
    to_process = []
    domain_to_uniprot = {}
    domain_to_uniprot_name = {}
    for pdb, annot in inp['annotations'].items():
        domain, domain_annot = next(iter(annot.items()))
        uniprot = domain_annot['uniprot_id']
        uniprot_name = domain_annot['uniprot_name']
        chain = domain_annot['chain_id']
        to_process.append((domain, pdb, chain, uniprot))
        domain_to_uniprot[domain] = uniprot
        domain_to_uniprot_name[domain] = uniprot_name
    
    # Download GO term information from PDBeAPI
    domains_to_terms = defaultdict(list)
    terms_to_description = {}
    for domain, pdb, chain, uniprot in to_process:
        response_str = requests.get(API_URL.format(pdb=pdb)).text
        response = json.loads(response_str)
        go_terms = response[pdb]['GO']
        for term, term_data in go_terms.items():
            for mapping in term_data['mappings']:
                if mapping['struct_asym_id'] == chain:
                    domains_to_terms[domain].append(term)
                    terms_to_description[term] = term_data['category'], term_data['name'], term_data['definition']
        print(domain, file=sys.stderr)

    result = {'terms': domains_to_terms, 'descriptions': terms_to_description, 'uniprot': domain_to_uniprot, 'uniprot_name': domain_to_uniprot_name}
    json.dump(result, sys.stdout)



if __name__ == '__main__':
    args = parse_args()
    exit_code = main(**args)
    if exit_code is not None:
        exit(exit_code)