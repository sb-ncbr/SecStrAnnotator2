# This Python3 script downloads the list of domains belonging to the specified Pfam family or CATH homologous superfamily and prints it in JSON in format { pdb: [[domain_name, chain, range]] }.

import sys
import requests
import json
import argparse

#  PARSE ARGUMENTS  ################################################################################

DEFAULT_API_URL = 'https://www.ebi.ac.uk/pdbe/api/mappings'  # 'https://wwwdev.ebi.ac.uk/pdbe/api/mappings'
DEFAULT_NUMBERING = 'label'

parser = argparse.ArgumentParser()
parser.add_argument('accession', help='"accession" argument to PDBe API SIFTS Mappings, e.g. Pfam accession or CATH cathcode.', type=str)
parser.add_argument('--numbering', help='Numbering scheme for residue numbers and chain IDs (default = ' +DEFAULT_NUMBERING + ').', choices=['label', 'auth'], default=DEFAULT_NUMBERING)
parser.add_argument('--source', help='URL with PDBeAPI server (default = ' + DEFAULT_API_URL + ')', default=DEFAULT_API_URL)
parser.add_argument('--allow_null_domain_name', help='Allow domain name to be null if not provided by the source (otherwise will create CATH-like domain names)', action='store_true')
parser.add_argument('--join_domains_in_chain', help='Join all domains in one chain if their names are not provided by the source', action='store_true')
parser.add_argument('--chain_change_warning', help='Print warning if struct_asym_id != chain_id (label_ vs auth_ chain numbering)', action='store_true')
args = parser.parse_args()

accession = args.accession
api_url = args.source
allow_null_domain_name = args.allow_null_domain_name
join_domains_in_chain = args.join_domains_in_chain
chain_change_warning = args.chain_change_warning


if args.numbering == 'label':
	AUTH_NUMBERING = False
elif args.numbering == 'auth':
	AUTH_NUMBERING = True
else:
	raise Exception('Unknown numbering: ' + args.numbering)

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

def get_start(mapping):
	return mapping['start']['author_residue_number'] if AUTH_NUMBERING else mapping['start']['residue_number']

def get_end(mapping):
	return mapping['end']['author_residue_number'] if AUTH_NUMBERING else mapping['end']['residue_number']

def get_chain(mapping):
	return mapping['chain_id'] if AUTH_NUMBERING else mapping['struct_asym_id']

def get_domains_multisegment(mappings, pdb):  # Returns a list of tuples (domain_name, chain, ranges).
	segments = []
	for mapping in mappings:
		domain = get_domain_name(mapping)
		chain = get_chain(mapping)
		start = get_start(mapping)
		end = get_end(mapping)
		if chain_change_warning and mapping['chain_id'] != mapping['struct_asym_id']:
			sys.stderr.write(f'Warning: {pdb} {mapping["chain_id"]} --> {mapping["struct_asym_id"]}\n')
		safe_str = lambda i: str(i) if i is not None else ''
		rang = safe_str(start) + ':' + safe_str(end)
		segments.append((domain, chain, rang))
	if all( domain is not None for domain, chain, rang in segments ):
		# All segments have a domain name, join segments with the same domain name.
		segments_by_domain = create_multidict( (domain, (chain, rang)) for domain, chain, rang in segments )
		result = []
		for domain, chain_ranges in sorted(segments_by_domain.items()):
			chains, ranges = zip(*chain_ranges)
			chains = set(chains)
			if len(chains) > 1:
				sys.stderr.write('Error: Some domains contain parts of multiple chains\n')
				exit(1)
			chain = next(iter(chains))
			ranges = ','.join( str(rang) for rang in ranges )
			result.append((domain, chain, ranges))
	else:
		# Some segments miss a domain name, creating new domain names.
		# sys.stderr.write('Some segments miss a domain name, creating new domain names\n')
		ranges_by_chain = create_multidict( (chain, rang) for domain, chain, rang in segments )
		result = []
		for chain, ranges in sorted(ranges_by_chain.items()):
			if allow_null_domain_name and join_domains_in_chain:
				result.append( (None, chain, ','.join(ranges)) )
			elif allow_null_domain_name:
				result.extend( (None, chain, rang) for rang in ranges )
			elif len(ranges) == 1 or join_domains_in_chain:
				result.append( (pdb + chain + '00', chain, ','.join(ranges)) )
			else:
				result.extend( (pdb + chain + str(i).zfill(2), chain, rang) for i, rang in enumerate(ranges, start=1) )
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
for pdb, entry in results.items():
	if isinstance(entry, list):
		output[pdb] = get_domains_multisegment(entry, pdb)
	elif isinstance(entry, dict) and 'mappings' in entry:
		output[pdb] = get_domains_multisegment(entry['mappings'], pdb)

n_pdbs = len(output)
n_domains = sum( len(doms) for pdb, doms in output.items() )

sys.stderr.write(f'Found {n_domains} domains in {n_pdbs} PDB entries.\n')

print(json.dumps(output, sort_keys=True, indent=4))