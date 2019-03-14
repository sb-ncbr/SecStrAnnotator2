# This Python3 script downloads macromolecular structures in PDB format.

import os
import sys
import json
import requests
import gzip
import argparse
import re
from collections import OrderedDict

#  CONSTANTS  ##############################################################################

URL_PDB = 'http://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pdb}.ent'
URL_PDB_GZ = 'http://ftp.ebi.ac.uk/pub/databases/rcsb/pdb-remediated/data/structures/divided/pdb/{pdb2}{pdb3}/pdb{pdb}.ent.gz'
URL_CIF = 'http://www.ebi.ac.uk/pdbe/entry-files/download/{pdb}_updated.cif'
URL_CIF_GZ = 'https://files.rcsb.org/download/{pdb}.cif.gz'

#  FUNCTIONS  ###############################################################################

def log(message=''):
	sys.stderr.write(message + '\n')

def try_read_json(filename):
	with open(filename) as f:
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
	with open(filename) as f:
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
		with open(output_file, 'w') as w:
			w.write(text)
		return True
	else: 
		return False

class ProgressBar:
	def __init__(self, n_steps, width=100, title='', writer=sys.stdout):
		self.n_steps = n_steps # expected number of steps
		self.width = width
		self.title = (' '+title+' ')[0:min(len(title)+2, width)]
		self.writer = writer
		self.done = 0 # number of completed steps
		self.shown = 0 # number of shown symbols
	def start(self):
		self.writer.write('|' + self.title + '_'*(self.width-len(self.title)) + '|\n')
		self.writer.write('|')
		self.writer.flush()
		return self
	def step(self, n_steps=1):
		if self.n_steps == 0:
			return
		self.done = min(self.done + n_steps, self.n_steps)
		new_shown = int(self.width * self.done / self.n_steps)
		self.writer.write('*' * (new_shown-self.shown))
		self.writer.flush()
		self.shown = new_shown
	def finalize(self):
		self.step(self.n_steps - self.done)
		self.writer.write('|\n')
		self.writer.flush()

#  PARSE ARGUMENTS  ###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='File with the list of PDBs to download', type=str)
parser.add_argument('output_directory', help='Directory to save downloaded files', type=str)
parser.add_argument('--input_format', 
	help=('Specify the input file format: '
		+ 'json = JSON file in format [PDB] or {PDB:anything}; '
		+ 'json_by_uniprot = JSON file in format {UniProtId:[PDB]} or {UniProtId:{PDB:anything}}; '
		+ 'text = text file with a list of PDBs (with any separator)'), 
	choices=['json', 'json_by_uniprot', 'text'],
	default='json')
parser.add_argument('--format', 
	help=('Specify the format of dowloaded structures'), 
	choices=['cif', 'pdb'],
	default='cif')
parser.add_argument('--unique_uniprot', help='Download only the first PDB entry for each UniProtID', action='store_true')
parser.add_argument('--no_gzip', help='Download uncompressed files instead of default gzip', action='store_true')
args = parser.parse_args()

if args.unique_uniprot and args.input_format != 'json_by_uniprot':
	sys.stderr.write('Error: "--unique_uniprot" can only be used with "--input_format json_by_uniprot" \n')

#  MAIN  ###############################################################################

# Check and prepare output directory
outdir = args.output_directory
if os.path.isfile(outdir):
	raise Exception('Cannot create output directory '+outdir+' because a file with this name exists.')
if not os.path.isdir(outdir):
	os.makedirs(outdir)

log('Reading file "' + args.input_file + '" in format "' + args.input_format + '"' + (', unique UniProtID.' if args.unique_uniprot else '.'))
log()

# Read input file
if args.input_format == 'json':
	pdbs = read_pdbs_json(args.input_file)
elif args.input_format == 'json_by_uniprot':
	pdbs = read_pdbs_json_by_uniprot(args.input_file, unique_uniprot=args.unique_uniprot)
elif args.input_format == 'text':
	pdbs = read_pdbs_text(args.input_file)

# Download files
failed = []
bar = ProgressBar(len(pdbs), title = 'Downloading ' + str(len(pdbs)) + ' PDBs').start()
for pdb in pdbs:
	ok = download_pdb(pdb, args.format, os.path.join(outdir, pdb + '.' + args.format), not args.no_gzip)
	if not ok:
		failed.append(pdb)
	bar.step()
bar.finalize()

n_failed = len(failed)
n_ok = len(pdbs) - n_failed
log()
log(f'Downloaded {n_ok} PDB entries, failed to download {n_failed} PDB entries')
if n_failed > 0:
	log('Failed to download: ' + ' '.join(failed))
	exit(1)

