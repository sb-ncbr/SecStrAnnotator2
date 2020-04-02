import os
import sys
import json
import requests
import argparse
import re
from collections import defaultdict
from typing import Tuple
from os import path

import lib

#  CONSTANTS  ##############################################################################

from constants import *

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('domain_list', help='JSON file with the list of domains in SecStrAPI format (from domain_lists_to_SecStrAPI_format.py)', type=str)
args = parser.parse_args()

domain_list_file = args.domain_list

#  MAIN  ##############################################################################

URL_TEMPLATE = 'https://www.uniprot.org/uniprot/{uniprot}.fasta'
RE_GENE=re.compile('GN=(\w+)')
RE_SPECIES=re.compile('OS=(\w+ \w+)')

with open(domain_list_file) as r:
    input_json = json.load(r)

pdb2domains = input_json[ANNOTATIONS]

uniprots = sorted(set(domain[UNIPROT_ID] for domains in pdb2domains.values() for domain in domains.values()))

for uniprot in uniprots:
    response = requests.get(URL_TEMPLATE.format(uniprot=uniprot)).text
    label, sequence = response.split('\n', maxsplit=1)
    if not label.startswith('>') or not sequence.replace('\n', '').isupper():
        raise Exception(f'Unexpected format:\n{response}')
    _, _, description = label[1:].split('|')
    gene_match = RE_GENE.search(description)
    gene = gene_match.group(1) if gene_match is not None else 'N/A'
    species_match = RE_SPECIES.search(description)
    species = species_match.group(1) if species_match is not None else 'N/A'
    print(f'>{uniprot} {species} {gene}\n{sequence}')

