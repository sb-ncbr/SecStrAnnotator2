# This file contains a PyMOL script that adds a new command 'annotate_sec_str'.
# To run this script once, type 'run FILENAME' into PyMOL (where FILENAME is the exact location of this file). The will also write out usage information.
# To run this script every time you launch PyMOL, add line 'run FILENAME' (where FILENAME is the exact location of this file) into file ~/.pymolrc (Linux, MacOS X) or C:\Program Files\DeLano Scientific\PyMOL\pymolrc (or C:\Program Files\PyMOL\PyMOL\pymolrc) (Windows).
# Version: 1.1 (2017/10/26)
# Further information about SecStrAPI is available at http://webchem.ncbr.muni.cz/API/SecStr/

# Constants
THIS_SCRIPT = 'annotate_sec_str_extension-1.0.py'
CHANNELSDB_API = 'https://webchem.ncbr.muni.cz/API/ChannelsDB/PDB/' # If no annotation file is provided, annotation will be downloaded from this URL (with PDB ID added after slash).
SEC_STR_API = 'http://webchem.ncbr.muni.cz/API/SecStr/Annotation/' # If no annotation file is provided, annotation will be downloaded from this URL (with PDB ID added after slash).
DEFAULT_BASE_COLOR = 'gray80' # This color will be used for residues without annotation.
DEFAULT_REPRESENTATION = 'cartoon' # Default visual representation for newly fetched and loaded structures (does not apply objects that have already been loaded).
DEFAULT_HET_REPRESENTATION = 'sticks' # Default visual representation for heteroatoms.
SHOW_LABELS = True # Indicates whether labels with SSE names will be shown.
LABEL_SIZE = None # Size for the labels (None = default size).

# Imports
import sys
import json
import requests
import re

from pymol import cmd

COLORS = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan', 'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine', 'olive', 'gray90', 'gray50', 'palegreen', 'lightblue', 'purpleblue']
global color_i
color_i = 0

# Field names in JSON annotation format
PDB = 'pdb'
SSES = 'secondary_structure_elements'
LABEL = 'label'
START_RESI = 'start'
END_RESI = 'end'
CHAIN_ID = 'chain_id'
TYPE = 'type'
SHEET_ID = 'sheet_id'
COLOR='color'

# Auxiliary functions

def test():
	fail('urmumma\nis fat')
	Annotation_1_0('')
	return

def annss(selection, annotation_file=None, name=None, base_color = DEFAULT_BASE_COLOR, force_cartoon=True):
	if annotation_file is not None and name is not None:
		annotation_text = get_annotation_from_file(annotation_file)
		annotation = parse_annotation_text(annotation_text)
		found_pdbids = annotation.get_pdbids()
		if name in found_pdbids:
			pdbid = name
			domains = annotation.get_domains_by_pdbid(pdbid)
		else:
			domains = annotation.get_domains_by_name(name)
			if len(domains) == 0:
				fail('Annotation for ' + name + ' not found in file ' + annotation_file)
				return False
			else:
				pdbid = domains[0].get_pdbid()
	elif name is not None:
		pdbid = name[0:4]
		annotation_text = get_annotation_from_internet(pdbid)
		annotation = parse_annotation_text(annotation_text)
		domains = annotation.get_domains_by_name(name)
	elif annotation_file is not None:
		pdbids = annotation.get_pdbids()
		if len(pdbids) == 1:
			pdbid = pdbids[0]
			log('Selected PDB ID: ' + pdbid + ' (the only one in the annotation file)')
		else:
			pdbid = first_occurring_sample(pdbids, selection)
			if pdbid is None:
				fail('Could not guess PDB ID from selection "' + selection + '" (PDB IDs in annotation file: ' + ', '.join(pdbids) + ')')
				return False
			else:
				log('Selected PDB ID: ' + pdbid)
		name = pdbid
		domains = annotation.get_domains_by_pdbid(pdbid)
	else: # guess pdbid using selection and download annotation
		words = re.split('\W+', selection)
		debug_log(words)
		try:
			name = next( word for word in words if is_valid_pdbid(word) or is_valid_domain_name(word) )
			debug_log(name)
			pdbid = first_valid_pdbid(name)
			log('Guessed PDB ID: ' + pdbid)
		except:
			fail('Could not guess PDB ID from selection "' + selection + '"')
			return False
		annotation_text = get_annotation_from_internet(pdbid)
		annotation = parse_annotation_text(annotation_text)
		domains = annotation.get_domains_by_pdbid(pdbid)

	log('PDB ID: ' + pdbid + ', NAME: ' + name)
	
	if not is_valid_selection(selection):
		fetch_structure(pdbid, force_cartoon=force_cartoon, domains=domains)
		if not is_valid_selection(selection):
			fail('Invalid selection "' + selection + '"')
			return False

	apply_annotation(selection, domains)
	return True


def debug_log(message):
	print(message)

def log(message):
	print(message)

def fail(message):
	sys.stderr.write(THIS_SCRIPT + ':\n    Error: ' + message.replace('\n', '\n    ') + '\n')
	# TODO print some helpfull info

def is_valid_pdbid(string):
	return len(string) == 4 and string.isalnum() and string[0].isdigit()
	
def is_valid_domain_name(string):
	return len(string) >= 4 and is_valid_pdbid(string[0:4])
	
def is_valid_selection(selection):
	try:
		cmd.count_atoms(selection)
		return True
	except:
		return False

def first_valid_pdbid(text):
	match = re.search('[0-9][a-zA-Z0-9]{3}', text)
	if match is not None:
		return match.group(0)
	else:
		return None

def first_occurring_sample(samples, text):
	indices_samples = [ (text.find(sample), sample) for sample in samples if sample in text ]
	if len(indices_samples) > 0:
		index, sample = min(indices_samples)
		return index
	else:
		return None

def fetch_structure(pdbid, force_cartoon=True, domains=None):
	if not is_valid_selection(pdbid):
		log('Fetching ' + pdbid)
		cmd.fetch(pdbid, async=0)
	if not is_valid_selection(pdbid):
		return False
	if force_cartoon:
		cmd.hide('everything', pdbid)
		if domains is None:
			sel = pdbid
		else:
			sel = ' or '.join( '(' + domain.get_selection() + ')' for domain in domains )
		cmd.show(DEFAULT_REPRESENTATION, sel)
		cmd.show(DEFAULT_HET_REPRESENTATION, '(' + sel + ') and het')
	return True

def get_annotation_from_file(filename):
	'''Reads annotation file in JSON format.'''
	with open(filename,'r') as f:
		annotation = f.read()
	return annotation

def get_annotation_from_internet(pdbid):
	pdbid = pdbid.lower()
	url = SEC_STR_API + pdbid
	log('Downloading annotation for ' + pdbid + '\n  [' + url + ']')
	try:
		annotation_text = requests.get(url).text
		debug_log(annotation_text)
		return annotation_text
	except:
		fail('Failed to download annotation for ' + pdbid)
		return None

def parse_annotation_text(text):
	return Annotation_1_0(text)

class Annotation_1_0:
	def __init__(self, text):
		try:
			self.json = json.loads(text)
		except:
			fail('Could not parse annotation')
		# TODO continue here
		raise NotImplementedError()
	def get_pdbids(self):
		raise NotImplementedError()
	def get_domains_by_pdbid(self, pdbid):
		raise NotImplementedError()
	def get_domains_by_name(self, name):
		raise NotImplementedError()

class Domain_1_0:
	def __init__(self, text):
		raise NotImplementedError()
	def get_pdbid(self):
		raise NotImplementedError()






# OLDER:

def annotate_sec_str(selection, annotation_file = None, base_color = BASE_COLOR):
	if not is_valid_selection(selection):
		pdb_id = selection[0:4]
		if is_valid_pdbid(pdb_id):
			if os.path.isfile(pdb_id+'.pdb'):
				print('Loading '+pdb_id+'.pdb')
				cmd.load(pdb_id+'.pdb')
			else:
				print('Fetching '+pdb_id)
				cmd.fetch(pdb_id,async=0)
			if not is_valid_selection(selection):
				print('\''+selection+'\' is not a valid selection')
				print('FAILED')
				return
			cmd.hide('everything', pdb_id)
			cmd.show(DEFAULT_REPRESENTATION, pdb_id)
			cmd.show(DEFAULT_HET_REPRESENTATION, pdb_id+' & het')
		else:
			print('\''+pdb_id+'\' is not a valid PDB ID')
			print('FAILED')	
			return
	object_names = cmd.get_names('objects',0,selection)
	for pdb_id in object_names:
		color_by_annotation(pdb_id, pdb_id+' & ('+selection+')', base_color, annotation_file)
	print('DONE')

def color_by_annotation(pdb_id, selection, base_color, annotation_file):
	'''Obtains annotation from an annotation file or from internet (if annotation_file==None) and use it to color the structure.'''
	if annotation_file!=None and annotation_file!='':
		sses = read_annotation_file_json(annotation_file, pdb_id)
	else:
		sses = download_annotation_json(pdb_id)
	cmd.color(base_color, selection+' & not het & symbol c')
	if SHOW_LABELS and LABEL_SIZE != None:
		cmd.set('label_size', str(LABEL_SIZE))
	for sse in sorted(sses, key = lambda x: x[CHAIN_ID]+str(x[START_RESI]).zfill(10)):
		label = sse[LABEL]
		chain_id = sse[CHAIN_ID]
		sel_name = pdb_id+chain_id+'_'+label.replace('\'','.') # apostrophe is not allowed in PyMOL identifiers
		sel_definition = '(' + selection + ') & chain '+sse[CHAIN_ID]+' & resi '+str(sse[START_RESI])+'-'+str(sse[END_RESI])+' & symbol c'
		if cmd.count_atoms(sel_definition)>0:
			cmd.select(sel_name, sel_definition)
			cmd.group('SSEs_'+pdb_id, sel_name, 'add')
			cmd.color(assign_color(sse), sel_name)
			if SHOW_LABELS:
				label_sse(sse, sel_name)
		cmd.deselect()

def read_annotation_file_json(filename, pdb_id):
	'''Reads annotation file in JSON format.'''
	with open(filename,'r') as f:
		annotation = json.load(f)
	if pdb_id in annotation:
		return annotation[pdb_id][SSES]
	else:
		print('\''+filename+'\' does not contain annotation for '+pdb_id)
		return []

def download_annotation_json(pdb_id):
	'''Downloads annotation from SecStrAPI.'''
	pdb_id = pdb_id.lower()
	print('Downloading annotation for '+pdb_id+' from '+SEC_STR_API+pdb_id)
	annotation = json.loads(requests.get(SEC_STR_API+pdb_id).text)
	if pdb_id in annotation:
		return annotation[pdb_id][SSES]
	else:
		print('No available annotation for '+pdb_id)
		return []

def assign_color(sse):
	sse_type = 'H' if (sse[TYPE] in 'GHIh') else 'E' if (sse[TYPE] in 'EBe') else ' '
	return 's'+str(hash(sse[LABEL])%1000).zfill(3)

def assign_color(sse):
	sse_type = 'H' if (sse[TYPE] in 'GHIh') else 'E' if (sse[TYPE] in 'EBe') else ' '
	if COLOR in sse:
		return sse[COLOR]
	else:
		return 's'+str(hash(sse[LABEL])%1000).zfill(3)

def label_sse(sse, selection):
	middle = (sse[START_RESI] + sse[END_RESI]) // 2
	label_selection = '(' + selection + ') & name cb & resi ' + str(middle)
	if cmd.count_atoms(label_selection)==0:
		label_selection = '(' + selection + ') & name ca & resi ' + str(middle)
	cmd.label(label_selection, '\''+sse[LABEL].replace('\'','\\\'')+'\'')

# The script for PyMOL:
# test()

# cmd.extend('annotate_sec_str', annotate_sec_str)
cmd.extend('annss', annss)
# print('')
# print('Command "annotate_sec_str" have been added.')
# print('')
# print('  Usage: annotate_sec_str selection [, annotation_file [, base_color]]')
# print('    default base_color: '+BASE_COLOR)
# print('    default annotation_file: None (downloaded from SecStrAPI instead)')
# print('  Examples: annotate_sec_str 1tqn')
# print('            annotate_sec_str 1og2 and chain a, my_annotation.json, white')
# print('')

