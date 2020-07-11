'''
This Python3 script downloads full amino acids sequences for listed domains and prints them as multi-FASTA.

Example usage:
    python3  get_full_sequences.py  annotations.json 
'''

import argparse
from typing import Dict, Any, Optional
import os
import sys
import json
import requests
import math

import lib
from constants import *

#  CONSTANTS  ################################################################################

URL_TEMPLATE = 'https://www.ebi.ac.uk/proteins/api/proteins/{uniprot}'
# URL_TEMPLATE = 'https://www.uniprot.org/uniprot/{uniprot}.fasta'

#  FUNCTIONS  ################################################################################

def get_gene_name(api_response):
    try:
        return api_response['gene'][0]['name']['value']
    except:
        pass
    try:
        return api_response['gene'][0]['olnNames'][0]['value']
    except:
        pass
    try: 
        return api_response['gene'][0]['orfNames'][0]['value']
    except:
        pass
    return ''

def get_species_name(api_response):
    species = api_response['organism']['names'][0]['value']
    species = ' '.join(species.split()[:2])
    return species

def wrap_sequence(sequence: str, width=60) -> str:
    n_lines = math.ceil(len(sequence) / width)
    return '\n'.join(sequence[width*i : width*(i+1)] for i in range(n_lines))

#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('domain_list_file', help='JSON file with the list of domains in SecStrAPI format (from domain_lists_to_SecStrAPI_format.py)', type=str)
    args = parser.parse_args()
    return vars(args)


def main(domain_list_file: str) -> Optional[int]:
    '''Download full amino acids sequences for listed domains and print them as multi-FASTA.'''

    with open(domain_list_file, 'r', encoding=lib.DEFAULT_ENCODING) as r:
        input_json = json.load(r)

    pdb2domains = input_json[ANNOTATIONS]

    uniprots = sorted(set(domain[UNIPROT_ID] for domains in pdb2domains.values() for domain in domains.values()))

    progress_bar = lib.ProgressBar(len(uniprots), title=f'Downloading sequences for {len(uniprots)} UniProt IDs', writer=sys.stderr).start()
    for uniprot in uniprots:
        url = URL_TEMPLATE.format(uniprot=uniprot)
        response = requests.get(url).text
        response = json.loads(response)
        if 'errorMessage' in response:
            sys.stderr.write(f'Response from {url} contains error: {response["errorMessage"]}\n')
        sequence = response['sequence']['sequence']
        gene = get_gene_name(response)
        species = get_species_name(response)
        print(f'>{uniprot} {species} {gene}\n{wrap_sequence(sequence)}\n')
        progress_bar.step()
    progress_bar.finalize()


if __name__ == '__main__':
    args = parse_args()
    exit_code = main(**args)
    if exit_code is not None:
        exit(exit_code)