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

# Field names in SecStrAPI format - general and domain specification
API_VERSION = 'api_version'
ANNOTATIONS = 'annotations'
PDB = 'pdb'
CHAIN = 'chain'
RANGES = 'ranges'
UNIPROT_ID = 'uniprot_id'
UNIPROT_NAME = 'uniprot_name'
DOMAIN_MAPPINGS = 'domain_mappings'
DOMAIN = 'domain'
SOURCE = 'source'
FAMILY_ID = 'family'

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
				pdbid = domains[0][PDB]
	elif annotation_file is None and name is not None:
		pdbid = name[0:4]
		annotation_text = get_annotation_from_internet(pdbid)
		annotation = parse_annotation_text(annotation_text)
		domains = annotation.get_domains_by_name(name)
	elif annotation_file is not None and name is None:
		annotation_text = get_annotation_from_file(annotation_file)
		annotation = parse_annotation_text(annotation_text)
		pdbids = annotation.get_pdbids()
		if len(pdbids) == 1:
			pdbid = pdbids[0]
			log('Selected PDB ID: ' + pdbid + ' (the only one in the annotation file)')
		else:
			pdbid = first_valid_pdbid(selection, pdbids=pdbids)
			if pdbid is None:
				fail('Could not guess PDB ID from selection "' + selection + '" (PDB IDs in annotation file: ' + ', '.join(pdbids) + ')')
				return False
			else:
				log('Selected PDB ID: ' + pdbid)
		name = pdbid
		domains = annotation.get_domains_by_pdbid(pdbid)
	else: # annotation_file and name are None, guess pdbid using selection and download annotation
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
		pdbid_as_in_selection = str(first_valid_pdbid(selection, pdbids=[pdbid]))  # in case that the PBDID is used in uppercase in the selection
		try:
			fetch_structure(pdbid_as_in_selection, force_cartoon=force_cartoon, domains=domains)
		except None: # TODO change to except:
			warn('Failed to fetch "' + pdbid_as_in_selection + '"')
		if not is_valid_selection(selection):
			fail('Invalid selection "' + selection + '"')
			return False

	debug_log('OK')
	# apply_annotation(selection, domains)
	return True


def debug_log(message):
	print(message)

def log(message):
	print(message)

def warn(message):
	sys.stderr.write(THIS_SCRIPT + ':\n    Warning: ' + message.replace('\n', '\n    ') + '\n')

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

def first_valid_pdbid(text, pdbids=None):
	if pdbids is None:
		match = re.search('[0-9][a-zA-Z0-9]{3}', text)
		if match is not None:
			return match.group(0)
		else:
			return None
	else:
		lower_text = text.lower()
		indices_lengths = [ (lower_text.find(pdbid.lower()), len(pdbid)) for pdbid in pdbids ]
		indices_lengths = [ (i, l) for (i, l) in indices_lengths if i >= 0 ]
		if len(indices_lengths) > 0:
			index, length = min(indices_lengths)
			return text[index:index+length]
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
			sel = ' or '.join( '(' + selection_from_domain(domain) + ')' for domain in domains )
			sel = pdbid + ' and (' + sel + ')'
		cmd.show(DEFAULT_REPRESENTATION, sel)
		cmd.show(DEFAULT_HET_REPRESENTATION, '(' + sel + ') and het')
	return True

def selection_from_domain(domain):
	selection = domain.get(PDB, 'all')
	if CHAIN in domain:
		selection += ' and chain ' + domain[CHAIN]
	if RANGES in domain:
		selection += ' and resi ' + ranges_to_pymol_style(domain[RANGES])
	return selection

def range_to_pymol_style(the_range):
	fro, to = the_range.split(':')
	return fro.replace('-', '\-') + '-' + to.replace('-', '\-')

def ranges_to_pymol_style(ranges):
	return '+'.join( range_to_pymol_style(rang) for rang in ranges.split(',') )

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
		debug_log(annotation_text[0:70] + '...')
		return annotation_text
	except:
		fail('Failed to download annotation for ' + pdbid)
		return None

def parse_annotation_text(text):
	try:
		js = json.loads(text)
	except:
		fail('Could not parse annotation (not a valid JSON)')
		raise
	version = js.get(API_VERSION, '0.9')
	if version == '0.9':
		return Annotation_0_9(js)
	elif version == '1.0':
		return Annotation_1_0(js)
	else:
		fail('The annotation file is of a newer version (' + version + ') than this script (1.0). Please download the updated version of this script from SecStrAPI website http://webchem.ncbr.muni.cz/API/SecStr.')
		raise NotImplementedError()

class Annotation_1_0:
	def __init__(self, json_object):
		self.pdb2domains = {}
		# Change domain organization for each PDBID from dict-like to list-like:
		for pdbid, domains in json_object.get(ANNOTATIONS, {}).items():
			if isinstance(domains, dict):
				for name, domain in domains.items():
					domain[DOMAIN] = name
				self.pdb2domains[pdbid] = list(domains.values())
			else:
				self.pdb2domains[pdbid] = list(domains)
	def get_pdbids(self):
		return sorted(self.pdb2domains.keys())
	def get_domains_by_pdbid(self, pdbid):
		return self.pdb2domains.get(pdbid.lower(), [])
	def get_domains_by_name(self, name):
		result = []
		for pdbid, domains in self.pdb2domains.items():
			for domain in domains:
				match = domain.get(DOMAIN, None) == name or any ( mapping[DOMAIN] == name for mapping in domain.get(DOMAIN_MAPPINGS, []) )
				if match:
					result.append(domain)
		return result

class Annotation_0_9:
	def __init__(self, json_object):
		self.domain_dict = json_object
		for name, domain in self.domain_dict.items():
			domain[DOMAIN] = name
			domain[PDB] = name[0:4]
	def get_pdbids(self):
		return sorted(set( domain[PDB] for domain in self.domain_dict.values() ))
	def get_domains_by_pdbid(self, pdbid):
		return [ domain for domain in self.domain_dict.values() if domain[PDB] == pdbid.lower() ]
	def get_domains_by_name(self, name):
		if name in self.domain_dict:
			return [self.domain_dict[name]]
		else:
			return []



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
print('')
print('Command "annss" has been added.')
print('')
print('annss(selection, annotation_file=None, name=None, base_color = DEFAULT_BASE_COLOR, force_cartoon=True)')
# TODO add usefull usage info
# print('  Usage: annotate_sec_str selection [, annotation_file [, base_color]]')
# print('    default base_color: '+BASE_COLOR)
# print('    default annotation_file: None (downloaded from SecStrAPI instead)')
# print('  Examples: annotate_sec_str 1tqn')
# print('            annotate_sec_str 1og2 and chain a, my_annotation.json, white')
# print('')

