# This Python3 script downloads the list of domains belonging to the specified Pfam family or CATH homologous superfamily and prints it in JSON in format { pdb: [[domain_name, chain, range]] }.

import sys
import requests
import json
import argparse

from constants import *

#  PARSE ARGUMENTS  ################################################################################

DEFAULT_API_URL = 'https://www.ebi.ac.uk/pdbe/api/mappings'  # 'https://wwwdev.ebi.ac.uk/pdbe/api/mappings'
DEFAULT_NUMBERING = 'label'

parser = argparse.ArgumentParser()
parser.add_argument('accession', help='"accession" argument to PDBe API SIFTS Mappings, e.g. Pfam accession or CATH cathcode.', type=str)
# parser.add_argument('--numbering', help='Numbering scheme for residue numbers and chain IDs (default = ' +DEFAULT_NUMBERING + ').', choices=['label', 'auth'], default=DEFAULT_NUMBERING)
parser.add_argument('--source', help='URL with PDBeAPI server (default = ' + DEFAULT_API_URL + ')', default=DEFAULT_API_URL)
# parser.add_argument('--allow_null_domain_name', help='Allow domain name to be null if not provided by the source (otherwise will create CATH-like domain names)', action='store_true')
parser.add_argument('--join_domains_in_chain', help='Join all domains in one chain if their names are not provided by the source', action='store_true')
parser.add_argument('--chain_change_warning', help='Print warning if struct_asym_id != chain_id (label_ vs auth_ chain numbering)', action='store_true')
args = parser.parse_args()

accession = args.accession
# numbering = args.numbering
api_url = args.source
# allow_null_domain_name = args.allow_null_domain_name
join_domains_in_chain = args.join_domains_in_chain
chain_change_warning = args.chain_change_warning


#  FUNCTIONS  ################################################################################

def create_multidict(key_value_pairs):  # Creates a dictionary with a list of multiple values per key. 
	multidict = {}
	for key, value in key_value_pairs:
		if not key in multidict:
			multidict[key] = []
		multidict[key].append(value)
	return multidict

def get_domain_name(mapping, default=None):
	return mapping.get('domain', default)

def get_start(mapping, numbering=DEFAULT_NUMBERING):
	return mapping['start']['author_residue_number'] if numbering == 'auth' else mapping['start']['residue_number']

def get_end(mapping, numbering=DEFAULT_NUMBERING):
	return mapping['end']['author_residue_number'] if numbering == 'auth' else mapping['end']['residue_number']

def get_chain(mapping, numbering=DEFAULT_NUMBERING):
	return mapping['chain_id'] if numbering == 'auth' else mapping['struct_asym_id']

def get_range(mapping, numbering=DEFAULT_NUMBERING):
	if numbering == 'auth':
		start = str(mapping['start']['author_residue_number'] or '') + mapping['start']['author_insertion_code']
		end = str(mapping['end']['author_residue_number'] or '') + mapping['end']['author_insertion_code']
	else:
		start = str(mapping['start']['residue_number'] or '')
		end = str(mapping['end']['residue_number'] or '')
	return start + ':' + end

def get_domains_multisegment(mappings, pdb):  # Returns a list of tuples (domain_name, chain, ranges).
	segments = []
	for mapping in mappings:
		domain_name = get_domain_name(mapping)
		chain = get_chain(mapping, 'label')
		rang = get_range(mapping, 'label')
		auth_chain = get_chain(mapping, 'auth')
		auth_rang = get_range(mapping, 'auth')
		if chain_change_warning and chain != auth_chain:
			sys.stderr.write(f'Warning: {pdb} {auth_chain} --> {chain}\n')
		segments.append({DOMAIN_NAME: domain_name, CHAIN: chain, RANGES: rang, AUTH_CHAIN: auth_chain, AUTH_RANGES: auth_rang})
	if all( segment[DOMAIN_NAME] is not None for segment in segments ):
		# All segments have a domain name, join segments with the same domain name.
		segments_by_domain = create_multidict( (seg[DOMAIN_NAME], seg) for seg in segments )
		result = []
		for domain_name, domain_segments in sorted(segments_by_domain.items()):
			chains = [seg[CHAIN] for seg in domain_segments]
			ranges = [seg[RANGES] for seg in domain_segments]
			auth_chains = [seg[AUTH_CHAIN] for seg in domain_segments]
			auth_ranges = [seg[AUTH_RANGES] for seg in domain_segments]
			chains = sorted(set(chains))
			auth_chains = sorted(set(auth_chains))
			if len(chains) > 1:
				sys.stderr.write('Error: Some domains contain parts of multiple chains\n')
				exit(1)
			chain = chains[0]
			auth_chain = auth_chains[0]
			ranges = ','.join( str(rang) for rang in ranges )
			auth_ranges = ','.join( str(rang) for rang in auth_ranges )
			result.append({DOMAIN_NAME: domain_name, CHAIN: chain, RANGES: ranges, AUTH_CHAIN: auth_chain, AUTH_RANGES: auth_ranges})
	else:
		# Some segments miss a domain name, creating new domain names.
		# sys.stderr.write('Some segments miss a domain name, creating new domain names\n')
		ranges_by_chain = create_multidict( (seg[CHAIN], seg) for seg in segments )
		result = []
		for chain, chain_segments in sorted(ranges_by_chain.items()):
			auth_chain = next(seg[AUTH_CHAIN] for seg in chain_segments)
			ranges = [seg[RANGES] for seg in chain_segments]
			auth_ranges = [seg[AUTH_RANGES] for seg in chain_segments]
			if join_domains_in_chain:
				ranges = ','.join( str(rang) for rang in ranges )
				auth_ranges = ','.join( str(rang) for rang in auth_ranges )
				result.append( {DOMAIN_NAME: None, CHAIN: chain, RANGES: ranges, AUTH_CHAIN: auth_chain, AUTH_RANGES: auth_ranges} )
			else:
				for rang, auth_rang in zip(ranges, auth_ranges):
					result.append( {DOMAIN_NAME: None, CHAIN: chain, RANGES: rang, AUTH_CHAIN: auth_chain, AUTH_RANGES: auth_rang} )
	return result

#  MAIN  ###############################################################################

url = api_url + '/' + accession
sys.stderr.write('Downloading ' + url + '\n')
response = requests.get(url)
if response.ok:
	results = json.loads(response.text).get(accession, {}).get('PDB', {})
else: 
	sys.stderr.write('HTTP request failed, status code ' + str(response.status_code) + '.\n')
	exit(1)

output = {}
for pdb, entry in sorted(results.items()):
	if isinstance(entry, list):
		output[pdb] = get_domains_multisegment(entry, pdb)
	elif isinstance(entry, dict) and 'mappings' in entry:
		output[pdb] = get_domains_multisegment(entry['mappings'], pdb)

n_pdbs = len(output)
n_domains = sum( len(doms) for pdb, doms in output.items() )

sys.stderr.write(f'Found {n_domains} domains in {n_pdbs} PDB entries.\n')

json.dump(output, sys.stdout, indent=4)