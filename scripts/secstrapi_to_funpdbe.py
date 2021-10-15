import json
import os
from os import path
import sys
import argparse
from datetime import datetime
from typing import Tuple, List, Dict
from collections import defaultdict

from funpdbe_validator.validator.validator import Validator
from funpdbe_validator.validator.residue_index import ResidueIndexes

################################################################################

USE_AUTH_NUMBERING = True
MERGE_SITES_BY_LABEL = True
INCLUDE_GENERIC_NUMBERING = False
DECIDE_CONFIDENCE_FROM_METRIC = True

SOFTWARE_VERSION = 'SecStrAnnotator 2.2'  # SecStrAnnotator version TODO put automatically assigned full version of SecStrAnnotator

DATA_RESOURCE = 'SecStrAPI'
RESOURCE_ENTRY_URL_BASE = 'https://webchem.ncbr.muni.cz/API/SecStr/Annotation/'
FIRST_SITE_ID = 1
REFERENCE_RESIDUE_GENERIC_NUMBER = 50

TODAY = datetime.today().strftime('%d/%m/%Y')

ECO = [{ 
    "eco_term": "computational combinatorial evidence used in automatic assertion",
    "eco_code": "ECO_0000246"
}]  # TODO what to put here? "computational combinatorial evidence used in automatic assertion" or "structural similarity evidence" --> question for Misi
CONFIDENCE_CLASSIFICATION = 'medium'  # TODO what to put here? medium/null? --> question for Misi

HIGH_CONFIDENCE_METRIC_THRESHOLD = 10.0
MEDIUM_CONFIDENCE_METRIC_THRESHOLD = 20.0

################################################################################

def create_funpdbe(pdb: str, pdb_annotation: dict, resource_version: str) -> dict:
    funpdbe = {}
    funpdbe['data_resource'] = DATA_RESOURCE
    funpdbe['resource_version'] = resource_version  # @optional
    funpdbe['software_version'] = SOFTWARE_VERSION  # @optional
    funpdbe['resource_entry_url'] = RESOURCE_ENTRY_URL_BASE + pdb  # @optional
    funpdbe['release_date'] = TODAY  # @optional
    funpdbe['pdb_id'] = pdb
    chains, sites = create_funpdbe_chains_and_sites(pdb, pdb_annotation)
    funpdbe['chains'] = chains
    funpdbe['sites'] = sites
    # funpdbe['additional_entry_annotations'] = additional_entry_annotations  # @optional
    funpdbe['evidence_code_ontology'] = ECO
    return funpdbe

def create_funpdbe_chains_and_sites(pdb: str, pdb_annotation: dict) -> Tuple[list, list]:
    if MERGE_SITES_BY_LABEL:
        sites = {}  # keyed by label
        chain2res2sites = defaultdict(lambda: defaultdict(lambda: []))  # e.g. {'A': {123: [2, 7]}}
        # create sites and link residues to sites
        for domain, domain_annotation in pdb_annotation.items():
            for sse in domain_annotation['secondary_structure_elements']:
                site_label = create_funpdbe_site_label(domain, sse)
                if site_label in sites:
                    # merge with existing site
                    site_id = sites[site_label]['site_id']
                else:
                    # create new site
                    site_id = len(sites) + FIRST_SITE_ID
                    new_site = create_funpdbe_site(site_id, pdb, domain, sse)
                    sites[site_label] = new_site
                chain = sse['chain_id']
                metric = sse['metric_value']
                for resi in range(sse['start'], sse['end']+1):
                    chain2res2sites[chain][resi].append((site_id, metric))
        sites = list(sites.values())
    else:
        sites = []
        site_id = FIRST_SITE_ID
        chain2res2sites = defaultdict(lambda: defaultdict(lambda: []))  # e.g. {'A': {123: [2, 7]}}
        # create sites and link residues to sites
        for domain, domain_annotation in pdb_annotation.items():
            for sse in domain_annotation['secondary_structure_elements']:
                new_site = create_funpdbe_site(site_id, pdb, domain, sse)
                chain = sse['chain_id']
                for resi in range(sse['start'], sse['end']+1):
                    chain2res2sites[chain][resi].append((site_id, metric))
                sites.append(new_site)
                site_id += 1
    # create chains
    chains = [ create_funpdbe_chain(chain_id, res2sites, sites) for chain_id, res2sites in chain2res2sites.items() ]
    return chains, sites

def create_funpdbe_site(site_id: int, pdb: str, domain: str, sse: dict) -> dict:
    site = {}
    site['site_id'] = site_id
    site['label'] = create_funpdbe_site_label(domain, sse)
    # site['source_database'] = "pdb"  # @optional
    # site['source_accession'] = pdb  # @optional
    # site['source_release_date'] = TODAY  # @optional
    # site['additional_site_annotations'] = additional_site_annotations  # @optional
    reference_resi = sse.get('pivot_residue', None)
    if INCLUDE_GENERIC_NUMBERING and reference_resi is not None:
        site['temp'] = {'reference_residue': reference_resi, 'simple_label': sse['label']}  # to be removed afterwards
    return site

def create_funpdbe_site_label(domain: str, sse: dict) -> str:
    label = sse['label']
    typ = sse['type']
    if typ in 'hGHI':
        return f'Helix {label}'
    elif typ in 'eEB':
        return f'Strand {label}'
    else:
        raise Exception(f'Unknown SSE type "{typ}"')

def create_funpdbe_chain(chain_id: str, res2sites: dict, sites: List[dict]) -> dict:
    chain = {}
    chain['chain_label'] = chain_id
    # chain['additional_chain_annotations'] = additional_chain_annotations  # @optional
    chain['residues'] = [ create_funpdbe_residue(resi, site_ids, sites) for resi, site_ids in res2sites.items() ]
    return chain

def create_funpdbe_residue(resi: int, site_ids: List[Tuple[int, float]], sites: List[dict]) -> dict:
    residue = {}
    residue['pdb_res_label'] = str(resi)
    residue['aa_type'] = ''
    # residue['additional_residue_annotations'] = additional_residue_annotations  # @optional
    # TODO put info about generic residue numbering into 'additional_residue_annotations', like 'C.50'?  --> try to put it here, ask Misi if OK
    residue['site_data'] = [ create_funpdbe_site_datum(site_id, metric) for site_id, metric in site_ids ]
    if INCLUDE_GENERIC_NUMBERING:
        generic_numbers = []
        for site_id in site_ids:
            site = sites[site_id - FIRST_SITE_ID]
            ref_resi = site.get('temp', {}).get('reference_residue')
            simple_label = site.get('temp', {}).get('simple_label')
            if ref_resi is not None and simple_label is not None:
                generic = f'{simple_label}.{REFERENCE_RESIDUE_GENERIC_NUMBER + resi - ref_resi}'
                generic_numbers.append(generic)
            if len(generic_numbers) > 1:
                raise Exception(f'Residue {resi} has more than one generic numberings.')
            elif len(generic_numbers) == 1:
                if 'additional_residue_annotations' not in residue:
                    residue['additional_residue_annotations'] = {}
                residue['additional_residue_annotations']['generic_residue_number'] = generic_numbers[0]
    return residue

def create_funpdbe_site_datum(site_id: int, metric_value: float) -> dict:
    site_datum = {}
    site_datum['site_id_ref'] = site_id
    # site_datum['raw_score'] = raw_score  # @optional
    # site_datum['confidence_score'] = confidence_score  # @optional
    site_datum['confidence_classification'] = confidence_based_on_metric(metric_value) if DECIDE_CONFIDENCE_FROM_METRIC else CONFIDENCE_CLASSIFICATION  # no idea ahat to put here
    # site_datum['aa_variant'] = aa_variant  # @optional
    return site_datum

def confidence_based_on_metric(metric_value: float) -> str:
    if metric_value <= HIGH_CONFIDENCE_METRIC_THRESHOLD:
        return 'high'
    elif metric_value <= MEDIUM_CONFIDENCE_METRIC_THRESHOLD:
        return 'medium'
    else:
        return 'low'

def add_aa_type_to_funpdbe(funpdbe: dict, label2auth_file: str, convert_label2auth=False) -> dict:
    index = {}
    with open(label2auth_file) as r:
        for line in r:
            if not line.lstrip().startswith('#'):
                chain_id, resi, auth_chain_id, auth_resi, auth_ins_code, label_comp_id = line.strip().split('\t')
                if auth_ins_code == '?':
                    auth_ins_code = ''
                else:
                    print(f'WARNING: Residue with insertion code {auth_resi}{auth_ins_code} (not tested)')
                index[(chain_id, resi)] = (auth_chain_id, auth_resi, auth_ins_code, label_comp_id)
    for chain in funpdbe['chains']:
        chain_id = chain['chain_label']
        if len(chain['residues']) < 1:
            raise Exception(f'Chain {chain_id} contains no residues')
        for residue in chain['residues']:
            resi = residue['pdb_res_label']
            auth_chain_id, auth_resi, auth_ins_code, aa_type = index[(chain_id, resi)] 
            if convert_label2auth:
                residue['pdb_res_label'] = auth_resi + auth_ins_code
            residue['aa_type'] = aa_type
        if convert_label2auth:
            chain['chain_label'] = auth_chain_id

def validate_funpdbe(pdb: str, filename: str) -> None:
    validator = Validator(DATA_RESOURCE) # Same as in the JSON
    validator.load_schema("funpdbe_validator/data/funpdbe_schema.json")
    validator.load_json(filename)
    if not validator.basic_checks():
        print(f'{pdb}: NOK: Basic checks failed')
    elif not validator.validate_against_schema():
        print(f'{pdb}: NOK: Validation against schema failed')
    else:
        # Passed data validations
        residue_indexes = ResidueIndexes(validator.json_data)
        if not residue_indexes.check_every_residue():
            # Passed the index validation
            print(f'{pdb}: NOK: Index validation failed')
        else:
            print(f'{pdb}: OK')

def remove_all_temp(obj: object) -> None:
    if isinstance(obj, dict):
        obj.pop('temp', None)
        for nested in obj.values():
            remove_all_temp(nested)
    elif isinstance(obj, (list, tuple, set)):
        for nested in obj:
            remove_all_temp(nested)

###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('input_annotation_file', type=str)
parser.add_argument('label2auth_directory', type=str)
parser.add_argument('output_directory', type=str)
args = parser.parse_args()

input_file = args.input_annotation_file
label2auth_dir = args.label2auth_directory
output_dir = args.output_directory

################################################################################

with open(input_file) as f:
    js = json.load(f)
secstrapi_version = js['api_version']
annotations = js['annotations']

os.makedirs(output_dir, exist_ok=True)

for pdb, pdb_annotation in annotations.items():
    funpdbe = create_funpdbe(pdb, pdb_annotation, secstrapi_version)
    add_aa_type_to_funpdbe(funpdbe, path.join(label2auth_dir, pdb + '.label2auth.tsv'), convert_label2auth=USE_AUTH_NUMBERING)
    remove_all_temp(funpdbe)
    output_file = path.join(output_dir, pdb + '.json')
    with open(output_file, 'w') as w:
        json.dump(funpdbe, w, indent=4)
    validate_funpdbe(pdb, output_file)


# Insertion codes:
# > adam ~/Workspace/C#/SecStrAnnot2_data/SecStrAPI/testing_20190520/structures $ grep -e '^\w\+\W\w\+\W\w\+\W\w\+\W\w\+'  *.label2auth.tsv
# 6bld.label2auth.tsv:A	360	A	361	A
# 6bld.label2auth.tsv:A	361	A	361	B
