# PymolScriptSession for SecStrAnnotator 1.0
# Creates and saves a session with superimposed template and query structures and with selected and colored SSEs.
# Usage: 
#     pymol -qcyr show_alignment.py  directory  template_pdbid,template_chain,template_ranges  query_pdbid,query_chain,query_ranges  alignment_file


# Modifiable constants

TEMPLATE_BASE_COLOR       = 'wheat'  # This color will be used for template residues without annotation.
QUERY_BASE_COLOR          = 'white'  # This color will be used for query residues without annotation.
PROTEIN_REPRESENTATIONS   = ['ribbon', 'labels']  # Default visual representations for polymer.
LIGAND_REPRESENTATIONS    = ['sticks', 'labels']  # Default visual representations for heteroatoms.
BOND_REPRESENTATIONS      = ['dashes']  # Default visual representations for hydrogen bonds.
SHOW_LABELS               = True  # Indicates whether labels with SSE names will be shown.
LABEL_SIZE                = 30  # Size for the labels (None = default size).
LABEL_OUTLINE_COLOR       = None  # Outline color for the labels (None = default).
SYMBOL_INSTEAD_APOSTROPHE = '.'  # Apostrophe is not allowed in PyMOL identifiers and will be changed for this symbol


# Imports

from os import path 
import sys
import json
from pymol import cmd


# Get command line arguments

N_PSEUDOARGUMENTS = 1
arguments = sys.argv[N_PSEUDOARGUMENTS:]
options = [arg for arg in arguments if arg.startswith('-')]
arguments = [arg for arg in arguments if not arg.startswith('-')]
if len(arguments) == 4:
	directory, template, query, alignment = arguments
else:
	print('ERROR: Exactly 4 command line arguments required, ' + str(len(arguments)) + ' given: ' + ' '.join(arguments))
	cmd.quit(1)


# Additional constants

SSES = 'secondary_structure_elements'
HBONDS = 'hydrogen_bonds'
LABEL = 'label'
START_RESI = 'start'
END_RESI = 'end'
CHAIN_ID = 'chain_id'
TYPE = 'type'
SHEET_ID = 'sheet_id'
COLOR='color'

TEMPLATE_STRUCT_EXT = '.cif'
QUERY_STRUCT_EXT = '.cif'
TEMPLATE_ANNOT_EXT = '-template.sses.json'
QUERY_ANNOT_EXT = '-annotated.sses.json'
SESSION_EXT = '.pse'
TEMPLATE_OBJECT_PREFIX = 'T_'
TEMPLATE_SELECTION_PREFIX = TEMPLATE_OBJECT_PREFIX
QUERY_OBJECT_PREFIX = ''
QUERY_SELECTION_PREFIX = QUERY_OBJECT_PREFIX
HBONDS_OBJECT_NAME = 'hbonds'


# Auxiliary functions

def read_annotation_file_json(filename, pdb_id):
	"""Reads annotation file in JSON format."""
	with open(filename,'r') as f:
		annotation = json.load(f)
	if pdb_id in annotation:
		return annotation[pdb_id]
	else:
		print('\''+filename+'\' does not contain annotation for '+pdb_id)
		return {}

def color_by_annotation(pdb_id, selection, base_color, annotation_file, selection_prefix=''):
	"""Obtains annotation from an annotation file and use it to color the structure."""
	annotation = read_annotation_file_json(annotation_file, pdb_id)
	sses = annotation.get(SSES, [])
	cmd.color(base_color, selection + ' & not het & symbol c')
	if SHOW_LABELS and LABEL_SIZE != None:
		cmd.set('label_size', str(LABEL_SIZE))
		cmd.set('label_color', base_color, selection)
		if LABEL_OUTLINE_COLOR is not None:
			cmd.set('label_outline_color', LABEL_OUTLINE_COLOR)
	for sse in sorted(sses, key = lambda x: (x[CHAIN_ID], x[START_RESI])):
		label = sse[LABEL]
		chain_id = sse[CHAIN_ID]
		start = sse[START_RESI]
		end = sse[END_RESI]
		sel_name = selection_prefix + safe_object_name(label)
		sel_expression = selection_expression(selection, chain_id, str(start) + ':' + str(end), symbol='C')
		if cmd.count_atoms(sel_expression) > 0:
			cmd.select(sel_name, sel_expression)
			cmd.color(assign_color(sse), sel_name)
			if SHOW_LABELS:
				label_sse(sse, selection)
	if show_hbonds:
		hbonds = annotation.get(HBONDS, [])
		distance_name = selection_prefix + HBONDS_OBJECT_NAME
		for donor_chain, donor_resi, acceptor_chain, acceptor_resi in hbonds:
			donor = '(%s) and chain %s and resi %s and name N' % (selection, donor_chain, donor_resi)
			acceptor = '(%s) and chain %s and resi %s and name O' % (selection, acceptor_chain, acceptor_resi)
			cmd.distance(distance_name, donor, acceptor)
	cmd.deselect()

def assign_color(sse):
	"""Assigns a color to, SSE based on 'color' field or generated from the label."""
	if COLOR in sse:
		return sse[COLOR]
	else:
		magic_value = hash(sse[LABEL])
		return 's' + str(magic_value % 1000).zfill(3)

def label_sse(sse, selection):
	"""Adds a label to SSE"""
	middle = (sse[START_RESI] + sse[END_RESI]) // 2
	label_selection = '(' + selection + ') & name cb & resi ' + str(middle)
	# print('Label selection: '+label_selection)
	if cmd.count_atoms(label_selection)==0:
		label_selection = '(' + selection + ') & name ca & resi ' + str(middle)
	cmd.label(label_selection, '\''+sse[LABEL].replace('\'','\\\'')+'\'')

def safe_object_name(name):
	"""Replaces apostrophes, which are not allowed by PyMOL."""
	return name.replace('\'', SYMBOL_INSTEAD_APOSTROPHE)  # apostrophe is not allowed in PyMOL identifiers

def convert_ranges(ranges):
	"""Converts changes from SecStrAnnotator format (1:10,15:20) to PyMOL format (1-10+15-20)."""
	return ranges.replace('-','\-').replace(':','-').replace(',','+')

def selection_expression(object_name, chain, ranges, symbol=None):
	"""Creates PyMOL selection expression."""
	parts = [object_name]
	if chain is not None:
		parts.append('chain ' + chain)
	if ranges is not None:
		parts.append('resi ' + convert_ranges(ranges))
	if symbol is not None:
		parts.append('symbol ' + symbol)
	return ' & '.join(parts) + ' '

def apply_representations(selection, protein_reprs=[], ligand_reprs=[], bonds_reprs=[]):
	"""Hides everything of selection and shows given representations for HET and NOT HET atoms."""
	protein_sel = '(' + selection + ') & not het'
	ligand_sel = '(' + selection + ') & het'
	cmd.hide('everything', selection)
	for repr in protein_reprs:
		cmd.show(repr, protein_sel)
	for repr in ligand_reprs:
		cmd.show(repr, ligand_sel)
	for repr in bonds_reprs:
		cmd.show(repr, '*'+HBONDS_OBJECT_NAME)

def pdb_chain_ranges(domain_spec):
	parts = domain_spec.split(',', 2)
	pdb = parts[0] 
	chain = parts[1] if len(parts) >= 2 else None
	ranges = parts[2] if len(parts) >= 3 else None
	return (pdb, chain, ranges)

def read_alignment(filename):
    with open(filename) as f:
        text = f.read()
    try:
        aln = json.loads(text)
        return [(chain1, str(resi1), chain2, str(resi2)) for chain1, resi1, chain2, resi2 in aln['aligned_residues']]
    except Exception:
        lines = [line for line in text.split('\n') if len(line) > 0 and not line.startswith('#')]
        aln = []
        for line in lines:
            chain1, resi1, _, _, _, resi2, chain2 = line.split('\t')  # A	10	K	0.757	P	15	A
            match = (chain1, resi1, chain2, resi2)
            if '-' not in match:
                aln.append(match)
        return aln


def create_session(directory, template, query, alignment_file):
    """Creates and saves a session with superimposed template and query and selected and colored SSEs."""
    files = []

    template_id, template_chain, template_range = pdb_chain_ranges(template)  # template.split(',', 2)
    template_struct_file = path.join(directory, template_id + TEMPLATE_STRUCT_EXT)
    files.append(template_struct_file)

    query_id, query_chain, query_range = pdb_chain_ranges(query)  # query.split(',', 2)
    query_struct_file = path.join(directory, query_id + QUERY_STRUCT_EXT)
    session_file = path.join(directory, 'aln-'+template_id+'-'+query_id + SESSION_EXT)
    files.append(query_struct_file)
	
    for filename in files:
        if not path.isfile(filename):
            print('ERROR: File not found: ' + filename)
            cmd.quit(1)

    template_object = TEMPLATE_OBJECT_PREFIX + template_id
    cmd.load(template_struct_file, template_object)
    
    query_object = QUERY_OBJECT_PREFIX + query_id
    cmd.load(query_struct_file, query_object)
    
    cmd.hide('everything', 'all')
    apply_representations(selection_expression(template_object, template_chain, None), protein_reprs=PROTEIN_REPRESENTATIONS, ligand_reprs=LIGAND_REPRESENTATIONS)
    apply_representations(selection_expression(query_object, query_chain, None), protein_reprs=PROTEIN_REPRESENTATIONS, ligand_reprs=LIGAND_REPRESENTATIONS)
    alignment = read_alignment(alignment_file)
    for chain1, resi1, chain2, resi2 in alignment:
        sel1 = '%s and chain %s and resi %s and name CA' % (template_object, chain1, resi1)
        sel2 = '%s and chain %s and resi %s and name CA' % (query_object, chain2, resi2)
        cmd.distance('dist', sel1, sel2)
    cmd.color('grey70', template_object)
    cmd.spectrum('resi', 'rainbow', query_object + ' and symbol C')
    cmd.hide('labels')
    cmd.dss()
    cmd.zoom('vis')
    cmd.save(session_file)


# Main script

try:
	cmd.set('cif_use_auth', False)

	create_session(directory, template, query, alignment)
	cmd.quit(0)
except Exception as e:
	print('ERROR: Exception raised: ' + str(e))
	cmd.quit(1)