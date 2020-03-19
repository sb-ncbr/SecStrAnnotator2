# This file contains a PyMOL script that adds a new command 'annotate_sec_str'.
# To run this script once, type 'run FILENAME' into PyMOL (where FILENAME is the exact location of this file). The will also write out usage information.
# To run this script every time you launch PyMOL, add line 'run FILENAME' (where FILENAME is the exact location of this file) into file ~/.pymolrc (Linux, MacOS X) or C:\Program Files\DeLano Scientific\PyMOL\pymolrc (or C:\Program Files\PyMOL\PyMOL\pymolrc) (Windows).
# Version: 1.1 (2017/10/26)
# Further information about SecStrAPI is available at http://webchem.ncbr.muni.cz/API/SecStr/

# Constants
CHANNELSDB_API = 'https://webchem.ncbr.muni.cz/API/ChannelsDB/PDB/' # If no annotation file is provided, annotation will be downloaded from this URL (with PDB ID added after slash).
SEC_STR_API = 'http://webchem.ncbr.muni.cz/API/SecStr/Annotation/' # If no annotation file is provided, annotation will be downloaded from this URL (with PDB ID added after slash).
BASE_COLOR = 'gray80' # This color will be used for residues without annotation.
DEFAULT_REPRESENTATION = 'cartoon' # Default visual representation for newly fetched and loaded structures (does not apply objects that have already been loaded).
DEFAULT_HET_REPRESENTATION = 'sticks' # Default visual representation for heteroatoms.
SHOW_LABELS = True # Indicates whether labels with SSE names will be shown.
LABEL_SIZE = None # Size for the labels (None = default size).

# Imports
import json
import requests

from pymol import cmd
from pymol.cgo import *
import chempy
def addAtom(model, name, vdw, x, y, z, partialCharge = 0.0):
	a = chempy.Atom()
	a.name = name
	a.vdw = vdw
	a.coord = [x, y, z]
	a.partial_charge = partialCharge
	model.atom.append(a)

CATEGORY_SHORTCUTS = {'CofactorTunnels': 'cofactor', 'CSATunnels': 'csa', 'ReviewedChannels': 'reviewed', 'TransmembranePores': 'pore'}
COLORS = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan', 'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine', 'olive', 'gray90', 'gray50', 'palegreen', 'lightblue', 'purpleblue']
global color_i
color_i = 0

def download_channels_json(pdb_id):
	'''Downloads channels from ChannelsDB.'''
	pdb_id = pdb_id.lower()
	print('Downloading channels for ' + pdb_id)
	print('[' + CHANNELSDB_API + pdb_id + ']')
	channel_json = json.loads(requests.get(CHANNELSDB_API + pdb_id).text)
	channels = channel_json.get('Channels', {})
	annotations = channel_json.get('Annotations', [])
	reviewed = channels.get('ReviewedChannels', [])
	for rev in reviewed:
		identifier = rev['Id']
		try:
			rev['Name'] = next(annot['Name'] for annot in annotations if annot['Id'] == identifier)
		except StopIteration:
			pass
	return channels

def show_channels(selection):
	pdb_id = selection[0:4]
	annotate_sec_str(selection)
	cmd.dss()
	cmd.hide('everything', 'not (' + selection + ')')
	channels = download_channels_json(pdb_id)
	group = 'Channels_' + pdb_id
	cmd.group(group)
	for category, channs in channels.items():
		if len(channs) > 0:
			subgroup = group + '.' + category
			cmd.group(subgroup)
		# short_cat = CATEGORY_SHORTCUTS.get(category, category)
		for channel in channs:
			show_channel(channel, group=subgroup)
		# cmd.group(group, category, 'add')
	# cmd.group(group, 'open')
	cmd.zoom('vis')

def show_channel(channel, group=None):
	name = channel['Name'].replace(' ', '_') + '_' + channel['Id'] if 'Name' in channel else channel['Id']
	name = group + '.' + name
	global color_i
	color = COLORS[color_i]
	color_i = (color_i+1) % len(COLORS)
	profile = channel.get('Profile', [])
	print('    ' + name)
	model = chempy.models.Indexed()
	for i, point in enumerate(profile):
		addAtom(model, str(i), point['Radius'], point['X'], point['Y'], point['Z'])
	for a in range(len(model.atom)-1):
		b = chempy.Bond()
		b.index = [a,a+1]
		model.bond.append(b)
	cmd.set('surface_mode', 1)
	cmd.set('sphere_mode', 5)
	cmd.set('mesh_mode', 1)
	cmd.load_model(model, name)
	cmd.hide('everything', name)
	cmd.set('sphere_color', color, name)
	cmd.show('spheres', name)
	# cmd.label(label_selection, name)
	# cmd.group(group, name, 'add')

# Field names in JSON annotation format
SSES = 'secondary_structure_elements'
LABEL = 'label'
START_RESI = 'start'
END_RESI = 'end'
CHAIN_ID = 'chain_id'
TYPE = 'type'
SHEET_ID = 'sheet_id'
COLOR='color'

# Auxiliary functions
def annotate_sec_str(selection, annotation_file = None, base_color = BASE_COLOR):
	if not is_valid_selection(selection):
		pdb_id = selection[0:4]
		if is_valid_pdb_id(pdb_id):
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

def is_valid_pdb_id(pdb_id):
	return len(pdb_id)==4 and pdb_id.isalnum() and pdb_id[0].isdigit()
	
def is_valid_selection(selection):
	try:
		cmd.count_atoms(selection)
		return True
	except:
		return False

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
cmd.extend('annotate_sec_str', annotate_sec_str)
cmd.extend('show_channels', show_channels)
print('')
print('Commands "show_channels" and "annotate_sec_str" have been added.')
print('')
print('  Usage: show_channels selection')
print('  Examples: show_channels 1tqn')
print('            show_channels 1og2 and chain A')
print('')
print('  Usage: annotate_sec_str selection [, annotation_file [, base_color]]')
print('    default base_color: '+BASE_COLOR)
print('    default annotation_file: None (downloaded from SecStrAPI instead)')
print('  Examples: annotate_sec_str 1tqn')
print('            annotate_sec_str 1og2 and chain a, my_annotation.json, white')
print('')

