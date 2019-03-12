# Example of usage:
# python3 ../../SecStrAnnot2/scripts/merge_domain_lists.py CATH 1.10.630.10 cyps_cath_20190312.json Pfam PF00067 cyps_pfam_20190312.json --api_version 1.0 > cyps_merged_20190312.json

import sys
import json
import requests
import argparse
from collections import defaultdict
from typing import Tuple

#  CONSTANTS  ##############################################################################

DOMAINS_IN_DICT = False

API_VERSION = 'api_version'
ANNOTATIONS = 'annotations'
PDB = 'pdb'
CHAIN = 'chain'
RANGES = 'ranges'
UNIPROT_ID = 'uniprot_id'
UNIPROT_NAME = 'uniprot_name'
MAPPINGS = 'domain_mappings'
DOMAIN_NAME = 'domain'
SOURCE = 'source'
FAMILY_ID = 'family'

#  FUNCTIONS  ##############################################################################

class UniProtManager:
    UNIPROT = 'UniProt'
    UNIPROT_NAME = 'name'
    MAPPINGS = 'mappings'
    CHAIN = 'struct_asym_id'
    def __init__(self, source_url='http://www.ebi.ac.uk/pdbe/api/mappings/uniprot/'):
        self.pdbchain2uniidname = {}
        self.source_url = source_url
    def get_uniprot_id_and_name(self, pdb, chain) -> Tuple[str, str]:  
        if (pdb, chain) not in self.pdbchain2uniidname:
            self._add_pdb(pdb)
            if (pdb, chain) not in self.pdbchain2uniidname:
                self.pdbchain2uniidname[(pdb, chain)] = (None, None)
        return self.pdbchain2uniidname[(pdb, chain)]
    def _add_pdb(self, pdb):
        url = f'{self.source_url}/{pdb}'
        r = json.loads(requests.get(url).text)
        unip = r[pdb][self.UNIPROT]
        for uni_id, uni_annot in unip.items():
            uni_name = uni_annot[self.UNIPROT_NAME]
            for mapping in uni_annot[self.MAPPINGS]:
                chain = mapping[self.CHAIN]
                if (pdb, chain) in self.pdbchain2uniidname and self.pdbchain2uniidname[(pdb, chain)] != (uni_id, uni_name):
                    old_id, old_name = self.pdbchain2uniidname[(pdb, chain)]
                    raise Exception(f'{pdb} chain {chain} has an ambiguous mapping to Uniprot: {old_id} vs {uni_id}')
                else:
                    self.pdbchain2uniidname[(pdb, chain)] = (uni_id, uni_name)

def new_mapping(domain_name, source, family, pdb, chain, ranges):
    mapping = { DOMAIN_NAME: domain_name, SOURCE: source, FAMILY_ID: family, PDB: pdb, CHAIN: chain, RANGES: ranges }
    return mapping

def create_domain_id(pdb, chain):
    return pdb + chain

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('lists', help='One or more triples of arguments SOURCE FAMILY FILE, where SOURCE is the name of a source (CATH, Pfam...), FAMILY is the identifier of the family in the source (1.10.630.10, PF00067...), and FILE is a JSON file with the list of domains from the source in format {PDB:[[domain, chain, ranges]]}', nargs='*')
parser.add_argument('--api_version', help='API version information to include in the output', type=str, default=None)
args = parser.parse_args()

if len(args.lists) % 3 != 0:
    parser.error('The number of LISTS arguments must be a multiple of 3')
else:
    sources = args.lists[0::3]
    families = args.lists[1::3]
    files = args.lists[2::3]
api_version = args.api_version
if api_version is None:
    sys.stderr.write('WARNING: "api_version" is set to null\n')

#  MAIN  ##############################################################################

annotations = {}
uniprot_manager = UniProtManager()

for source, family, filename in zip(sources, families, files):
    with open(filename) as r:
        pdb2domains = json.load(r)
    for pdb, domains in pdb2domains.items():
        if pdb not in annotations:
            annotations[pdb] = {} if DOMAINS_IN_DICT else []
        for domain_name, chain, ranges in domains:
            if DOMAINS_IN_DICT:
                to_unify = [ dom for dom in annotations[pdb].values() if dom[CHAIN] == chain ]
            else:
                to_unify = [ dom for dom in annotations[pdb] if dom[CHAIN] == chain ]
            if len(to_unify) == 0:  # create a new domain
                our_domain_id = create_domain_id(pdb, chain)
                dom = { PDB: pdb, CHAIN: chain, RANGES: ':' }  # Take whole chain always
                dom[UNIPROT_ID], dom[UNIPROT_NAME] = uniprot_manager.get_uniprot_id_and_name(pdb, chain)
                # dom[TEMPLATE] = template # TODO implement in annotation script, not here
                dom[MAPPINGS] = [new_mapping(domain_name, source, family, pdb, chain, ranges)]
                if DOMAINS_IN_DICT:
                    annotations[pdb][our_domain_id] = dom
                else:
                    annotations[pdb].append(dom)
            elif len(to_unify) == 1:  # just add a mapping to an existing domain
                dom = to_unify[0]
                dom[MAPPINGS].append(new_mapping(domain_name, source, family, pdb, chain, ranges))
            else:
                raise Exception(f'Domain {domain_name} ({pdb} {chain} {ranges}) can be unified with more than one other domain')

annotations = { pdb: annot for pdb, annot in sorted(annotations.items()) }

result = { API_VERSION: api_version, ANNOTATIONS: annotations }
json.dump(result, sys.stdout, indent=4)
print()
