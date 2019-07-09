import os
import sys
import json
import requests
import argparse

#  CONSTANTS  ##############################################################################

from constants import *

URL = 'http://www.ebi.ac.uk/pdbe/api/pdb/entry/publications/'

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('domain_list', help='JSON file with the list of domains in SecStrAPI format (from merge_domain_lists.py)', type=str)
args = parser.parse_args()

domain_list_file = args.domain_list

#  MAIN  ##############################################################################

with open(domain_list_file, 'r') as r:
	input_annotations = json.load(r)[ANNOTATIONS]

publinfo = {}

for pdb in input_annotations:
	response = json.loads(requests.get(URL + pdb).text)
	if pdb in response:
		publinfo[pdb] = response[pdb]
	else:
		sys.stderr.write(f'Error: {pdb}\n')

print(json.dumps(publinfo, indent=4))