# This is a PyMOL extension which adds a new command 'annotate_sec_str'.
# To run once, type 'run THIS_FILE' into PyMOL (where THIS_FILE is the exact location of this file).
# To run every time you launch PyMOL, install via Plugin Manager in PyMOL or add line 'run THIS_FILE' into your pymolrc file.
# Further information will be available after installation by typing 'help annotate_sec_str'.
# Information about SecStrAPI is available at http://webchem.ncbr.muni.cz/API/SecStr/
# Version: 1.2 (2019/09/18)


# Settings ####################################################################
# Modify this section to change the default behaviour:

SEC_STR_API_URL              = 'http://webchem.ncbr.muni.cz/API/SecStr/Annotation/'  # If no annotation file is provided, annotation will be downloaded from this URL (with PDB ID added after slash).
DEFAULT_BASE_COLOR           = 'gray80'  # This color will be used for residues without annotation, unless specified by parameter base_color.
DEFAULT_REPRESENTATION       = 'cartoon'  # Default visual representation for newly fetched structures (does not apply objects that have already been loaded).
DEFAULT_HET_REPRESENTATION   = 'sticks'  # Default visual representation for heteroatoms.
SHOW_LABELS                  = True  # Indicates whether labels with SSE names should be shown.
LABEL_SIZE                   = None  # Size for the labels (None = default size).
PYMOL_REPLACEMENT_CHARACTER  = '+'  # In selection names, avoid-characters will be replaced by PYMOL_REPLACEMENT_CHARACTER (i.e. all except alphanumeric characters A-Z a-z 0-9 and characters _ . + -).


# Imports #####################################################################

import sys
import json
import requests
import re
from os import path
from collections import OrderedDict
from pymol import cmd


# Constants ###################################################################

THIS_SCRIPT = 'annotate_sec_str_extension-1.2.py'

# Field names in SecStrAPI format 
API_VERSION = 'api_version'
ANNOTATIONS = 'annotations'
PDB = 'pdb'
CHAIN = 'chain_id'
AUTH_CHAIN = 'auth_chain_id'
RANGES = 'ranges'
AUTH_RANGES = 'auth_ranges'
UNIPROT_ID = 'uniprot_id'
UNIPROT_NAME = 'uniprot_name'
DOMAIN_MAPPINGS = 'domain_mappings'
DOMAIN = 'domain'
SOURCE = 'source'
FAMILY_ID = 'family'
SSES = 'secondary_structure_elements'
ROTATION_MATRIX = 'rotation_matrix'

# Per-SSE field names 
LABEL = 'label'
START_RESI = 'start'
AUTH_START_RESI = 'auth_start'
AUTH_START_RESI_INS_CODE = 'auth_start_ins_code'
END_RESI = 'end'
AUTH_END_RESI = 'auth_end'
AUTH_END_RESI_INS_CODE = 'auth_end_ins_code'
PIVOT_RESI = 'pivot_residue'
AUTH_PIVOT_RESI = 'auth_pivot_residue'
AUTH_PIVOT_RESI_INS_CODE = 'auth_pivot_residue_ins_code'
CHAIN_ID = 'chain_id'
AUTH_CHAIN_ID = 'auth_chain_id'
TYPE = 'type'
SHEET_ID = 'sheet_id'
COLOR='color'
INSERTION_CODE_NULL = '?'

# Auxiliary functions #########################################################

def debug_log(message):
	print(message)

def log(message):
	print(message)

def warn(message):
	sys.stderr.write(THIS_SCRIPT + ':\n')
	# sys.stderr.write('    Warning: ' + message.replace('\n', '\n    ') + '\n')
	print('    Warning: ' + message.replace('\n', '\n    '))

def fail(message):
	sys.stderr.write(THIS_SCRIPT + ':\n')
	# sys.stderr.write('    Error: ' + message.replace('\n', '\n    ') + '\n')
	print('    Error: ' + message.replace('\n', '\n    '))
	# TODO print some helpfull info

def is_valid_pdbid(string):
	return len(string) == 4 and string.isalnum() and string[0].isdigit()
	
def is_valid_domain_name(string):
	return len(string) >= 4 and is_valid_pdbid(string[0:4])
	
def is_valid_selection(selection):
	try:
		# debug_log('is_valid_selection(' + selection + ') ...')
		cmd.count_atoms(selection)
		# debug_log('is_valid_selection(' + selection + ') Yes')
		return True
	except:
		# debug_log('is_valid_selection(' + selection + ') No')
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

def fetch_structure(pdbid, domain_name=None, domain=None, force_cartoon=True):
	# debug_log('fetch structure ' + str(force_cartoon))
	if not is_valid_selection(pdbid):
		log('Fetching ' + pdbid)
		cmd.fetch(pdbid, async=0)
		if not is_valid_selection(pdbid):
			return False
		if force_cartoon:
			cmd.hide('everything', pdbid)
			if domain is None:
				sel = pdbid
			else:
				use_auth = should_use_auth(pdbid)
				# sel = ' or '.join( '(' + selection_from_domain(domain, use_auth) + ')' for domain in domains )
				sel = pdbid + ' and (' + selection_from_domain(domain, use_auth) + ')'
				if domain_name is not None and not is_valid_selection(domain_name):
					cmd.select(domain_name, sel)
			cmd.show(DEFAULT_REPRESENTATION, sel)
			cmd.show(DEFAULT_HET_REPRESENTATION, '(' + sel + ') and hetatm')
			cmd.zoom(sel)
	return True

def has_domain_auth_fields(domain):
	return all( domain.get(field, None) is not None for field in (AUTH_CHAIN, AUTH_RANGES) )

def selection_from_domain(domain, use_auth):
	if use_auth: 
		chain = domain.get(AUTH_CHAIN, None)
		ranges = domain.get(AUTH_RANGES, None)
	else:
		chain = domain.get(CHAIN, None)
		ranges = domain.get(RANGES, None)
	if chain is not None and ranges is not None:
		return '(chain ' + chain + ' and resi ' + ranges_to_pymol_style(ranges) + ')'
	elif chain is not None:
		return '(chain ' + chain + ')'
	else:
		return '(all)'

def has_sse_auth_fields(sse):
	return all( sse.get(field, None) is not None for field in (AUTH_CHAIN_ID, AUTH_START_RESI, AUTH_END_RESI) )

def selection_from_sse(sse, use_auth):
	if use_auth and any( sse.get(field, None) is None for field in (AUTH_CHAIN_ID, AUTH_START_RESI, AUTH_END_RESI) ):
		warn('auth_* fields are missing in the annotation for SSE "' + sse.get(LABEL, 'null') + '", using label_* fields instead')
		use_auth = False
	if use_auth:
		chain = sse[AUTH_CHAIN_ID]
		start = str(sse[AUTH_START_RESI]).replace('-', '\-')
		start_ins = sse.get(AUTH_START_RESI_INS_CODE, INSERTION_CODE_NULL)
		if start_ins != INSERTION_CODE_NULL:
			start += start_ins
		end = str(sse[AUTH_END_RESI]).replace('-', '\-')
		end_ins = sse.get(AUTH_END_RESI_INS_CODE, INSERTION_CODE_NULL)
		if end_ins != INSERTION_CODE_NULL:
			end += end_ins
	else:
		chain = sse[CHAIN_ID]
		start = str(sse[START_RESI]).replace('-', '\-')
		end = str(sse[END_RESI]).replace('-', '\-')
	selection = '(chain ' + chain + ' and resi ' + start + '-' + end + ')'
	return selection  
	# TODO solve insertion codes (resi 64-65A selects 64, 65, 65A, 65B, 65C!)
	# pepseq will not work when mapping annotation onto different structure
	# Or just say it's not my problem but problem of PyMOL.

def range_to_pymol_style(the_range):
	fro, to = the_range.split(':')
	return fro.replace('-', '\-') + '-' + to.replace('-', '\-')

def ranges_to_pymol_style(ranges):
	return '+'.join( range_to_pymol_style(rang) for rang in ranges.split(',') )

def is_from_cif(object_name):
	try:
		n_cif_atoms = cmd.count_atoms('(' + object_name + ') and not segi ""')
		return n_cif_atoms > 0
	except:
		return False

def should_use_auth(object_name):
	try:
		cif_use_auth = parse_boolean(cmd.get('cif_use_auth'))
	except:
		cif_use_auth = True
	return cif_use_auth or not is_from_cif(object_name)

def pymol_safe_name(name):
	chars = set(name)
	for char in chars:
		if not (char.isalnum() or char in '_.+-'):
			name = name.replace(char, PYMOL_REPLACEMENT_CHARACTER)
	return name

def get_annotation_from_file(filename):
	'''Reads annotation file in JSON format.'''
	with open(filename,'r') as f:
		annotation = f.read()
	return annotation

def get_annotation_from_internet(pdbid):
	pdbid = pdbid.lower()
	url = SEC_STR_API_URL + pdbid
	log('Downloading annotation for PDB entry ' + pdbid + '\n  [' + url + ']')
	try:
		annotation_text = requests.get(url).text
		# debug_log(annotation_text[0:70] + '...')
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
		fail('The annotation file is of a newer version (' + version + ') than this script supports (1.0). Please download the updated version of this script from SecStrAPI website http://webchem.ncbr.muni.cz/API/SecStr.')
		raise NotImplementedError()

class Annotation_1_0:
	def __init__(self, json_object):
		if 'error' in json_object:
			warn('Annotation file reports an error: ' + json_object['error'])
			self.pdb2domains = {}
			return
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
				match = name == domain.get(DOMAIN, None) or any ( name == mapping[DOMAIN] for mapping in domain.get(DOMAIN_MAPPINGS, []) )
				if match:
					result.append(domain)
		return result

class Annotation_0_9:
	def __init__(self, json_object):
		if 'error' in json_object:
			warn('Annotation file reports an error: ' + json_object['error'])
			self.domain_dict = {}
			return
		self.domain_dict = json_object
		for name, domain in self.domain_dict.items():
			domain[DOMAIN] = name
			domain[PDB] = name[0:4]
	def get_pdbids(self):
		return sorted(set( domain[PDB] for domain in self.domain_dict.values() ))
	def get_domains_by_pdbid(self, pdbid):
		return [ domain for domain in self.domain_dict.values() if pdbid.lower() == domain[PDB] ]
	def get_domains_by_name(self, name):
		if name in self.domain_dict:
			return [self.domain_dict[name]]
		else:
			return []

def apply_annotation(selection, domains, base_color, apply_rotation=False):
	# then = datetime.now()
	for domain in domains:
		selection = '(' + selection + ')'
		use_auth = should_use_auth(selection)

		if DOMAIN in domain:
			domain_name = domain[DOMAIN]
		else:
			domain_name = domain[PDB] + domain[CHAIN]
		domain_selection = selection + ' and ' + selection_from_domain(domain, use_auth)
		# debug_log('domain_name:"' + domain_name + '", domain selection: "' + domain_selection + '"')
		if not is_valid_selection(domain_name):  # such object may already be defined (e.g. if domain_name == pdb)
			cmd.select(domain_name, domain_selection)
		group = domain_name + '.sses'
		if not is_valid_selection(group):
			cmd.group(group)

		if apply_rotation and ROTATION_MATRIX in domain:
			print('Applying rotation for domain ' + domain_name)  # debug
			pdb = domain[PDB]
			matrix = domain[ROTATION_MATRIX]
			cmd.set_object_ttt(pdb, matrix)
			cmd.zoom(domain_name)

		if base_color is not None:
			cmd.color(base_color, domain_selection + ' and not hetatm and symbol C')
		if SHOW_LABELS and LABEL_SIZE != None:
			cmd.set('label_size', str(LABEL_SIZE))
		
		# debug_log('apply_annotation1: ' + str(datetime.now() - then))

		sses = domain.get(SSES, [])
		sses = sorted(sses, key=lambda sse: (sse[CHAIN_ID], sse[START_RESI]))
		if use_auth and not all( has_sse_auth_fields(sse) for sse in sses ):
			bad = [ sse.get(LABEL, 'null') for sse in sses if not has_sse_auth_fields(sse) ]
			bad = ('SSE ' + bad[0]) if len(bad) == 1 else ('SSEs ' + ', '.join(bad)) if len(bad) <= 5 else ('SSEs ' + ', '.join(bad[:5]) + '...')
			warn('auth_* fields are missing in the annotation (' + bad + '). Using label_* fields for domain "' + domain_name + '".')
			use_auth = False
		for sse in sses:
			label = sse[LABEL]
			safe_label = pymol_safe_name(label)
			chain_id = sse[CHAIN_ID]
			sel_name = group + '.' + safe_label
			sel_definition = selection + ' and ' + selection_from_sse(sse, use_auth)
			# debug_log('sse:"' + sel_name + '", sse selection: "' + sel_definition + '"')
			if True: # cmd.count_atoms(sel_definition) > 0: #debug
				cmd.select(sel_name, sel_definition)
				cmd.color(assign_color(sse), sel_name + ' and symbol C')
				if SHOW_LABELS:
					label_sse(sse, sel_name, use_auth)
			mark_pivot(sse, sel_name, use_auth, color=assign_color(sse))
			cmd.deselect()
	# debug_log('apply_annotation: ' + str(datetime.now() - then))

def assign_color(sse):
	sse_type = 'H' if (sse[TYPE] in 'GHIh') else 'E' if (sse[TYPE] in 'EBe') else ' '
	if COLOR in sse:
		return sse[COLOR]
	else:
		magic_number = hash(sse[LABEL]) % 1000
		return 's' + str(magic_number).zfill(3)

def label_sse(sse, selection, use_auth):
	if use_auth:
		chain = sse[AUTH_CHAIN_ID]
		start = sse[AUTH_START_RESI]
		end = sse[AUTH_END_RESI]
	else:
		chain = sse[CHAIN_ID]
		start = sse[START_RESI]
		end = sse[END_RESI]
	middle = (start + end) // 2
	label_selection = '(' + selection + ') and chain ' + chain + ' and resi ' + str(middle) + ' and name CB'
	if cmd.count_atoms(label_selection)==0:
		label_selection = '(' + selection + ') and chain ' + chain + ' and resi ' + str(middle) + ' and name CA'
	cmd.label(label_selection, '"' + sse[LABEL].replace('"', '\\"') + '"')

def mark_pivot(sse, selection, use_auth, color='red'):
	if use_auth:
		chain = sse[AUTH_CHAIN_ID]
		pivot = sse.get(AUTH_PIVOT_RESI, None)
		pivot_ins_code = sse.get(AUTH_PIVOT_RESI_INS_CODE, INSERTION_CODE_NULL)
		if pivot is not None:
			pivot = str(pivot)
			if pivot_ins_code != INSERTION_CODE_NULL:
				pivot += pivot_ins_code
	else:
		chain = sse[CHAIN_ID]
		pivot = sse.get(PIVOT_RESI, None)
		if pivot is not None:
			pivot = str(pivot)
	if pivot is not None:
		pivot_selection = '(' + selection + ') and chain ' + chain + ' and resi ' + pivot + ' and name CA+CB'
		cmd.show('sticks', pivot_selection)
		cmd.color(color, pivot_selection)
		pivot_selection = '(' + selection + ') and chain ' + chain + ' and resi ' + pivot + ' and name CA'
		cmd.show('spheres', pivot_selection)
		cmd.color('red', pivot_selection)

def parse_boolean(string):
	if not isinstance(string, str):
		return string
	if string.lower() in ['0', 'off', 'false', 'no']:
		return False
	elif string.lower() in ['1', 'on', 'true', 'yes']:
		return True
	else:
		fail('Invalid value "' + string + '", allowed values: 0, 1, off, on, false, true')
		raise Exception

# def print_help():
# 	print('annotate_sec_str(selection, annotation_file=None, name=None, base_color = DEFAULT_BASE_COLOR, force_cartoon=True)')
# 	print('Usage: annotate_sec_str selection [, annotation_file [, name [, base_color [, force_cartoon ]]]]')
# 	print('    default annotation_file: None (downloaded from SecStrAPI ' + SEC_STR_API_URL + ')')
# 	print('    default name: None (will be guessed from selection)')
# 	print('    default base_color: ' + DEFAULT_BASE_COLOR)
# 	print('    default force_cartoon: True')
# 	print('Examples: annotate_sec_str 1tqn')
# 	print('          annotate_sec_str 1og2, my_annotation.sses.json, 1og2A00')
# 	# TODO add usefull usage info

def test():
	base_dir = '/home/adam/Workspace/C#/SecStrAnnot2_data/SecStrAPI/testing_20190322'
	def run_test(expect_success, *args, **kwargs):
		log('---------------------------------------------------')
		log('TEST annotate_sec_str ' + ', '.join(str(arg) for arg in args) + ', ' + ', '.join(k + '=' + v for k, v in kwargs.items()) )
		cmd.delete('all')
		if annotate_sec_str(*args, **kwargs) == expect_success:
			log('TEST OK')
		else:
			log('TEST FAILED')
			raise Exception()
	log('---------------------------------------------------')
	# annot.file, name
	run_test(True, '1bu7', path.join(base_dir, 'annotations_with_pivots_all.json'), '1bu7')
	run_test(False, '1bu7', path.join(base_dir, 'annotations_with_pivots_best.json'), '1bu7')
	run_test(True, '1bu7', path.join(base_dir, 'annotations_with_pivots_all.json'), '1bu7A00')
	run_test(False, '1bu7', path.join(base_dir, 'annotations_with_pivots_all.json'), '1bu7xxx')
	# no annot.file, name
	run_test(True, '1bu7', None, '1bu7')
	# run_test(True, '1bu7', None, '1bu7A00')  # uncomment when annotations 1.0 are online
	run_test(False, '1bu7', None, '1bu7xxx')
	# annot.file, no name
	run_test(True, '1bu7', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(False, '1bu7_my_favourite', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(True, '1bu7A00', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(False, '1bu7A00_my_favourite', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(True, 'chain A and 1bu7', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(True, 'chain A and 1bu7A00', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(False, 'xxxx', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(False, '1xxx', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	# no annot.file, no name
	run_test(True, '1bu7', None, None)
	run_test(False, '1bu7_my_favourite', None, None)
	# run_test(True, '1bu7A00', None, None)  # uncomment when annotations 1.0 are online
	run_test(True, 'chain A and 1bu7', None, None)
	run_test(False, 'xxxx', None, None)
	run_test(False, '1xxx', None, None)

	# annot.file, name
	run_test(True, '1BU7', path.join(base_dir, 'annotations_with_pivots_all.json'), '1bu7')
	run_test(False, '1BU7', path.join(base_dir, 'annotations_with_pivots_best.json'), '1bu7')
	run_test(True, '1BU7', path.join(base_dir, 'annotations_with_pivots_all.json'), '1bu7A00')
	run_test(False, '1BU7', path.join(base_dir, 'annotations_with_pivots_all.json'), '1bu7xxx')
	# no annot.file, name
	run_test(True, '1BU7', None, '1bu7')
	# run_test(True, '1bu7', None, '1bu7A00')  # uncomment when annotations 1.0 are online
	run_test(False, '1BU7', None, '1bu7xxx')
	# annot.file, no name
	run_test(True, '1BU7', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(False, '1BU7_my_favourite', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(False, '1BU7A00', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(False, '1BU7A00_my_favourite', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(True, 'chain A and 1BU7', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(False, 'chain A and 1BU7A00', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(False, 'XXXX', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	run_test(False, '1XXX', path.join(base_dir, 'annotations_with_pivots_all.json'), None)
	# no annot.file, no name
	run_test(True, '1BU7', None, None)
	run_test(False, '1BU7_my_favourite', None, None)
	# run_test(True, '1bu7A00', None, None)  # uncomment when annotations 1.0 are online
	run_test(True, 'chain A and 1BU7', None, None)
	run_test(False, 'XXXX', None, None)
	run_test(False, '1XXX', None, None)

	cmd.delete('all')
	log('ALL TESTS OK')
	return True


# Main function ###############################################################

def annotate_sec_str(selection, annotation_file=None, name=None, base_color = DEFAULT_BASE_COLOR, force_cartoon=True):
	'''
DESCRIPTION

	"annotate_sec_str" visualizes protein secondary structure annotation.

USAGE

	annotate_sec_str selection, [, annotation_file [, name [, base_color [, force_cartoon ]]]]

ARGUMENTS

	selection = string: selection-expression. If it is not a valid selection-expression, the command will try to fetch the missing structures.

	annotation_file = string: name of the file with secondary structure annotation (in SecStrAPI format). If not specified, the annotation will be downloaded from SecStrAPI (the PDB ID will be guessed from the "name" argument).

	name = string: PDB ID or domain name specifying a protein domain(s) within the annotation file. If not specified, the PDB ID and the domain name will be guessed from the "selection" argument.

	base_color = string: color to use for non-annotated atoms {default: gray80}

	force_cartoon = 0/1: apply default cartoon representation for all newly fetched structures {default: 1}

NOTES

EXAMPLES

	annotate_sec_str 1tqn
	annotate_sec_str 1og2, my_annotation.sses.json, 1og2A00

	'''

	try:
		try:
			force_cartoon = parse_boolean(force_cartoon)
		except:
			return False
		fetched = False
		force_rotation = True
		# debug_log(force_cartoon)
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
					fail('Annotation for "' + name + '" not found')
					return False
				else:
					pdbid = domains[0][PDB]
		elif annotation_file is None and name is not None:
			pdbid = name[0:4]
			annotation_text = get_annotation_from_internet(pdbid)
			annotation = parse_annotation_text(annotation_text)
			domains = annotation.get_domains_by_name(name)
			if len(domains) == 0:
				fail('Annotation for "' + name + '" not found')
				return False
		elif annotation_file is not None and name is None:
			annotation_text = get_annotation_from_file(annotation_file)
			annotation = parse_annotation_text(annotation_text)
			words = re.split('\W+', selection)
			# debug_log(words)
			try:
				# name = next( word for word in words if len(annotation.get_domains_by_name(word, allow_prefix=True)) > 0 )
				name = next( word for word in words if len(annotation.get_domains_by_name(word)) > 0 )
				pdbid = name[0:4]
				# domains = annotation.get_domains_by_name(name, allow_prefix=True)
				domains = annotation.get_domains_by_name(name)
				if len(domains) == 0:
					fail('Annotation for "' + name + '" not found')
					return False
			except StopIteration:
				try:
					# name = next( word for word in words if len(annotation.get_domains_by_pdbid(word, allow_prefix=True)) > 0 )
					name = next( word for word in words if len(annotation.get_domains_by_pdbid(word)) > 0 )
					pdbid = name[0:4]
					# domains = annotation.get_domains_by_pdbid(name, allow_prefix=True)
					domains = annotation.get_domains_by_pdbid(name)
				except StopIteration:
					pdbids = annotation.get_pdbids()
					if len(pdbids) == 1:
						pdbid = pdbids[0]
						name = pdbid
						domains = annotation.get_domains_by_pdbid(pdbid)
						warn('Could not guess annotation name from selection "' + selection + '"')
						log('Selected PDB ID: "' + pdbid + '" (the only one in the annotation file)')
					else:
						fail('Could not guess annotation name from selection "' + selection + '"')
						return False
			# debug_log('Guessed PDB ID: ' + pdbid)
		else: # annotation_file and name are None, guess pdbid using selection and download annotation
			words = re.split('\W+', selection)
			# debug_log(words)
			try:
				name = next( word for word in words if is_valid_pdbid(word) or is_valid_domain_name(word) )
			except StopIteration:
				fail('Could not guess PDB ID from selection "' + selection + '"')
				return False
			pdbid = name[0:4]
			# log('Guessed PDB ID: ' + pdbid)
			annotation_text = get_annotation_from_internet(pdbid)
			annotation = parse_annotation_text(annotation_text)
			if name == pdbid:
				domains = annotation.get_domains_by_pdbid(pdbid)
			else:
				domains = annotation.get_domains_by_name(name)
			if len(domains) == 0:
				fail('Annotation for "' + name + '" not found')
				return False
		
		# TODO when guessing name, don't take next, but require unambiguity

		log('PDB ID: ' + pdbid + ', NAME: ' + name + ', DOMAINS: ' + str(len(domains)))
		# debug_log('Partial OK')
		
		if not is_valid_selection(selection):
			words = re.split('\W+', selection)
			for word in words:
				if is_valid_pdbid(word) and not is_valid_selection(word):
					fetch_structure(word, force_cartoon=force_cartoon)
					fetched = True
				elif is_valid_domain_name(word) and not is_valid_selection(word):
					pdb = word[0:4]
					if annotation_file is None:
						annotation_text = get_annotation_from_internet(pdb)
						annotation = parse_annotation_text(annotation_text)
					doms = annotation.get_domains_by_name(word)
					if len(doms) == 1:
						fetch_structure(pdb, domain_name=word, domain=doms[0], force_cartoon=force_cartoon)
						fetched = True
			if not is_valid_selection(selection):
				fail('Invalid selection "' + selection + '"')
				return False
		
		debug_log('OK')
		apply_annotation(selection, domains, base_color, apply_rotation=fetched and force_rotation)
		return True
	# except Exception as ex:
	# 	fail(str(ex))
	# 	return False
	except None:
		pass


# The script for PyMOL ########################################################

# test()
cmd.extend('annotate_sec_str', annotate_sec_str)
print('Command "annotate_sec_str" has been added. Type "help annotate_sec_str" for more information.')

#TODO rename pivot to reference, make mark_pivot optional


