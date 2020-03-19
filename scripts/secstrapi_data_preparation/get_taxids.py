import os
import sys
import json
import requests
import argparse

#  CONSTANTS  ##############################################################################

from constants import *

URL = 'http://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/'

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('domain_list', help='JSON file with the list of domains in SecStrAPI format (from merge_domain_lists.py)', type=str)
args = parser.parse_args()

domain_list_file = args.domain_list

#  MAIN  ##############################################################################

with open(domain_list_file, 'r') as r:
	input_annotations = json.load(r)[ANNOTATIONS]

for pdb, pdb_annot in input_annotations.items():
	response = json.loads(requests.get(URL + pdb).text)
	for domain, domain_annot in pdb_annot.items():
		chain = domain_annot[CHAIN]
		entities = [ entity for entity in response[pdb] if chain in entity['in_struct_asyms'] ] 
		taxids = [ source['tax_id'] for entity in entities for source in entity.get('source',[]) ]
		if len(set(taxids))!=1:
			sys.stderr.write('WARNING: ' + pdb + ' has non-unique taxon ID: ' + ', '.join([str(t) for t in taxids]) + '\n')
		print(domain, taxids[0], sep='\t')

exit()
