'''
This Python3 script downloads macromolecular structures in CIF or PDB format.

Example usage:
    python3  download_from_pdbe.py  domains.simple.json  my_structures/  --input_format json  --structure_format cif  --no_gzip  --cache cached_structures/
'''

import argparse
from typing import Dict, Any, Optional
import os
from os import path
import sys
import shutil
import json
import requests
import gzip
import re
from collections import OrderedDict

import lib

#  CONSTANTS  ################################################################################

URL_PDB = 'http://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pdb}.ent'
URL_PDB_GZ = 'http://ftp.ebi.ac.uk/pub/databases/rcsb/pdb-remediated/data/structures/divided/pdb/{pdb2}{pdb3}/pdb{pdb}.ent.gz'
URL_CIF = 'http://www.ebi.ac.uk/pdbe/entry-files/download/{pdb}_updated.cif'
URL_CIF_GZ = 'https://files.rcsb.org/download/{pdb}.cif.gz'
DEFAULT_INPUT_FORMAT = 'json'
DEFAULT_STRUCTURE_FORMAT = 'cif'

#  FUNCTIONS  ################################################################################

def log(message=''):
	sys.stderr.write(message + '\n')

def try_read_json(filename):
	with open(filename, 'r', encoding=lib.DEFAULT_ENCODING) as f:
		try:
			result = json.load(f, object_pairs_hook=OrderedDict)
		except ValueError as e:
			sys.stderr.write('Error: "' + filename + '" is not a valid JSON file (' + str(e) +  ') \n')
			exit(1)
	return result

def check_json_type(json_object, typeex):
	if isinstance(typeex, tuple):
		return any( check_json_type(json_object, typ) for typ in typeex )
	elif isinstance(typeex, str):
		return isinstance(json_object, str)
	elif isinstance(typeex, bool):
		return isinstance(json_object, bool)
	elif isinstance(json_object, (int, float)):
		return isinstance(json_object, (int, float))
	elif typeex is None:
		return json_object is None
	elif isinstance(typeex, list):
		if not isinstance(json_object, list):
			return False 
		if len(typeex) > 0:
			return all( check_json_type(elem, typeex[0]) for elem in json_object )
		else: 
			return True
	elif isinstance(typeex, dict):
		if not isinstance(json_object, dict):
			return False 
		if len(typeex) > 0:
			value_typeex = next(iter(typeex.values()))
			return all( isinstance(key, str) and check_json_type(value, value_typeex) for key, value in json_object.items() )
		else: 
			return True
	raise Exception('Invalid typeex')

def read_pdbs_text(filename):
	result = []
	with open(filename, 'r', encoding=lib.DEFAULT_ENCODING) as f:
		for line in iter(f.readline, ''):
			result.extend(word for word in re.split('\W+', line) if word != '')
	return result

def read_pdbs_json(filename):
	result = try_read_json(filename)
	if check_json_type(result, (['PDB'], {})):
		return list(result)
	else:
		sys.stderr.write('Error: Expected a JSON array or object in "' + filename + '" \n')
		exit(1)

def read_pdbs_json_by_uniprot(filename, unique_uniprot=False):
	pdbs_by_uniprot = try_read_json(filename)
	if check_json_type(pdbs_by_uniprot, ({'Uni':['PDB']}, {'Uni':{}})):
		result = []
		for uni, pdbs in pdbs_by_uniprot.items():
			if not isinstance(pdbs, (list, dict)):
				raise
			if unique_uniprot:
				first = next(iter(pdbs), None)
				if first is not None:
					result.append(first)
			else:
				result.extend(pdbs)
		return result
	else:
		sys.stderr.write('Error: Expected a JSON object containing arrays or objects in "' + filename + '" \n')
		exit(1)

def construct_url(pdb, file_format, use_gzip):
	if file_format == 'cif' and use_gzip:
		url = URL_CIF_GZ
	elif file_format == 'cif':
		url = URL_CIF
	elif file_format == 'pdb' and use_gzip:
		url = URL_PDB_GZ
	elif file_format == 'pdb':
		url = URL_PDB
	return url.format(pdb=pdb, pdb1=pdb[0], pdb2=pdb[1], pdb3=pdb[2], pdb4=pdb[3])

def download_pdb(pdb, file_format, output_file, use_gzip):
	url = construct_url(pdb, file_format, use_gzip)
	response = requests.get(url)
	if response.ok:
		if use_gzip:
			text = gzip.decompress(response.content).decode("utf-8")
		else:
			text = response.text
		with open(output_file, 'w', encoding=lib.DEFAULT_ENCODING) as w:
			w.write(text)
		return True
	else: 
		return False

def try_copy_pdb(pdb, file_format, output_file, cache_dir):
	if cache_dir is None:
		return False
	input_file = path.join(cache_dir, pdb + '.' + file_format)
	try:
		shutil.copyfile(input_file, output_file)
		return True
	except:
		return False

#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_file', help='File with the list of PDBs to download', type=str)
    parser.add_argument('output_directory', help='Directory to save downloaded files', type=str)
    parser.add_argument('--input_format', 
        help=(f'Specify the input file format (default: {DEFAULT_INPUT_FORMAT}): '
            + 'json = JSON file in format [PDB] or {PDB:anything}; '
            + 'json_by_uniprot = JSON file in format {UniProtId:[PDB]} or {UniProtId:{PDB:anything}}; '
            + 'text = text file with a list of PDBs (with any separator)'), 
        choices=['json', 'json_by_uniprot', 'text'],
        default=DEFAULT_INPUT_FORMAT)
    parser.add_argument('--structure_format', 
        help=(f'Specify the format of dowloaded structures (default: {DEFAULT_STRUCTURE_FORMAT})'), 
        choices=['cif', 'pdb'],
        default=DEFAULT_STRUCTURE_FORMAT)
    parser.add_argument('--unique_uniprot', help='Download only the first PDB entry for each UniProtID', action='store_true')
    parser.add_argument('--no_gzip', help='Download uncompressed files instead of default gzip', action='store_true')
    parser.add_argument('--cache', help='Directory from which the structures will by preferentially copied (if not found, falls back to download)', type=str)
    args = parser.parse_args()
    return vars(args)


def main(input_file: str, output_directory: str, input_format: str = DEFAULT_INPUT_FORMAT, structure_format: str = DEFAULT_STRUCTURE_FORMAT, 
        unique_uniprot: bool = False, no_gzip: bool = False, cache: Optional[str] = None) -> Optional[int]:
    '''Download macromolecular structures in CIF or PDB format.'''

    # Check and prepare output directory
    if os.path.isfile(output_directory):
        raise Exception(f'Cannot create output directory {output_directory} because a file with this name exists.')
    os.makedirs(output_directory, exist_ok=True)

    log(f'Reading file "{input_file}" in format "{input_format}"' + (', unique UniProtID' if unique_uniprot else '') + '.\n')

    # Read input file
    if input_format == 'json':
        pdbs = read_pdbs_json(input_file)
    elif input_format == 'json_by_uniprot':
        pdbs = read_pdbs_json_by_uniprot(input_file, unique_uniprot=unique_uniprot)
    elif input_format == 'text':
        pdbs = read_pdbs_text(input_file)

    # Download files
    failed = []
    bar = lib.ProgressBar(len(pdbs), title = f'Downloading {len(pdbs)} PDB structures', writer=sys.stderr).start()
    for pdb in pdbs:
        outfile = os.path.join(output_directory, pdb + '.' + structure_format)
        ok = try_copy_pdb(pdb, structure_format, outfile, cache)
        if not ok:
            ok = download_pdb(pdb, structure_format, outfile, not no_gzip)
        if not ok:
            failed.append(pdb)
        bar.step()
    bar.finalize()

    n_failed = len(failed)
    n_ok = len(pdbs) - n_failed
    log(f'\nDownloaded {n_ok} PDB entries, failed to download {n_failed} PDB entries')
    if n_failed > 0:
        log('Failed to download: ' + ' '.join(failed))
        return 1


if __name__ == '__main__':
    args = parse_args()
    exit_code = main(**args)
    if exit_code is not None:
        exit(exit_code)
