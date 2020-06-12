import os
import sys
import json
import requests
import argparse
import re
from collections import defaultdict
from typing import Tuple
from os import path
import math

import lib

#  CONSTANTS  ##############################################################################

from constants import *

#  FUNCTIONS  ###############################################################################

def get_gene_name(api_response):
    try:
        return response['gene'][0]['name']['value']
    except:
        pass
    try:
        return response['gene'][0]['olnNames'][0]['value']
    except:
        pass
    try: 
        return response['gene'][0]['orfNames'][0]['value']
    except:
        pass
    return ''

def get_species_name(api_response):
    species = response['organism']['names'][0]['value']
    species = ' '.join(species.split()[:2])
    return species

def wrap_sequence(sequence: str, width=60) -> str:
    n_lines = math.ceil(len(sequence) / width)
    return '\n'.join(sequence[width*i : width*(i+1)] for i in range(n_lines))

	
#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('domain_list', help='JSON file with the list of domains in SecStrAPI format (from domain_lists_to_SecStrAPI_format.py)', type=str)
args = parser.parse_args()

domain_list_file = args.domain_list

#  MAIN  ##############################################################################

# URL_TEMPLATE = 'https://www.uniprot.org/uniprot/{uniprot}.fasta'
URL_TEMPLATE = 'https://www.ebi.ac.uk/proteins/api/proteins/{uniprot}'

with open(domain_list_file) as r:
    input_json = json.load(r)

pdb2domains = input_json[ANNOTATIONS]

uniprots = sorted(set(domain[UNIPROT_ID] for domains in pdb2domains.values() for domain in domains.values()))

for uniprot in uniprots:
    response = requests.get(URL_TEMPLATE.format(uniprot=uniprot)).text
    response = json.loads(response)
    sequence = response['sequence']['sequence']
    gene = get_gene_name(response)
    species = get_species_name(response)
    print(f'>{uniprot} {species} {gene}\n{wrap_sequence(sequence)}\n')

