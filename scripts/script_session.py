# PymolScriptSession for SecStrAnnotator 1.0
# Creates and saves a session with superimposed template and query structures and with selected and colored SSEs.
# Usage: 
#     pymol -qcyr script_session.py -- file_format directory [ --hbonds ] [ template_pdbid,template_chain,template_ranges ] query_pdbid,query_chain,query_ranges


# Modifiable constants

# TEMPLATE_BASE_COLOR       = 'wheat'  # This color will be used for template residues without annotation.  # Optimal for session with black background
# QUERY_BASE_COLOR          = 'white'  # This color will be used for query residues without annotation.
TEMPLATE_BASE_COLOR       = 'brown'  # This color will be used for template residues without annotation.  # Better for image with white background
QUERY_BASE_COLOR          = 'gray50'  # This color will be used for query residues without annotation.

PROTEIN_REPRESENTATIONS   = ['cartoon', 'labels']  # Default visual representations for polymer.
LIGAND_REPRESENTATIONS    = ['sticks', 'labels']  # Default visual representations for heteroatoms.
LIGAND_EXTENSION          = 8  # Show ligands within this radius around the protein.
BOND_REPRESENTATIONS      = ['dashes']  # Default visual representations for hydrogen bonds.
SHOW_LABELS               = True  # Indicates whether labels with SSE names will be shown.
LABEL_SIZE                = 30  # Size for the labels (None = default size).
LABEL_OUTLINE_COLOR       = None  # Outline color for the labels (None = default).
SYMBOL_INSTEAD_APOSTROPHE = '+'  # Apostrophe is not allowed in PyMOL identifiers and will be changed for this symbol
IMAGE_SIZE                = (400, 400)  # Size of created PNG image (width, height) or None if no image should be saved


# Imports

from os import path 
import sys
import json
import random
from pymol import cmd


# Get command line arguments

N_PSEUDOARGUMENTS = 1
arguments = sys.argv[N_PSEUDOARGUMENTS:]
options = [arg for arg in arguments if arg.startswith('-')]
arguments = [arg for arg in arguments if not arg.startswith('-')]
if len(arguments) == 4:
    file_format, directory, template, query = arguments
elif len(arguments) == 3:
    file_format, directory, query = arguments
    template = None
else:
    print('ERROR: Exactly 3 or 4 command line arguments required, ' + str(len(arguments)) + ' given: ' + ' '.join(arguments))
    cmd.quit(1)
detected = '--detected' in options
show_hbonds = '--hbonds' in options


# Additional constants

if file_format == 'cif':
    USE_CIF = True
elif file_format == 'pdb':
    USE_CIF = False
else:
    raise Exception('Unknown file format: ' + file_format)

SSES = 'secondary_structure_elements'
HBONDS = 'hydrogen_bonds'
LABEL = 'label'
START_RESI = 'start'
END_RESI = 'end'
CHAIN_ID = 'chain_id'
TYPE = 'type'
SHEET_ID = 'sheet_id'
COLOR='color'

TEMPLATE_STRUCT_EXT = '.cif' if USE_CIF else '.pdb'
QUERY_STRUCT_EXT = ('-aligned' if template is not None else '') +  ('.cif' if USE_CIF else '.pdb')
TEMPLATE_ANNOT_EXT = '-template.sses.json'
QUERY_ANNOT_EXT = '-detected.sses.json' if detected else '-annotated.sses.json'
SESSION_EXT = '-detected.pse' if detected else '-annotated.pse'
TEMPLATE_PNG_EXT = '-template.png'
QUERY_PNG_EXT = '-detected.png' if detected else '-annotated.png'
TEMPLATE_OBJECT_PREFIX = ''
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
    elif len(annotation) == 1:
        key = next(iter(annotation.keys()))
        print('\''+filename+'\' does not contain annotation for '+pdb_id+', taking '+key+' instead')
        return annotation[key]
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
    cmd.group(selection_prefix)
    center_of_mass = cmd.centerofmass(selection)
    for sse in sorted(sses, key = lambda x: (x[CHAIN_ID], x[START_RESI])):
        label = sse[LABEL]
        if label.startswith('_'):
            continue
        chain_id = sse[CHAIN_ID]
        start = sse[START_RESI]
        end = sse[END_RESI]
        sel_name = selection_prefix + '.' + safe_object_name(label)
        sel_expression = selection_expression(selection, chain_id, str(start) + ':' + str(end), symbol='C')
        if cmd.count_atoms(sel_expression) > 0:
            cmd.select(sel_name, sel_expression)
            cmd.color(assign_color(sse), sel_name)
            if SHOW_LABELS:
                label_sse(sse, selection, center_of_mass=center_of_mass)
    if show_hbonds:
        hbonds = annotation.get(HBONDS, [])
        distance_name = selection_prefix + '.' + HBONDS_OBJECT_NAME
        for donor_chain, donor_resi, acceptor_chain, acceptor_resi in hbonds:
            donor = selection_expression(selection, donor_chain, donor_resi, name='N')
            acceptor = selection_expression(selection, acceptor_chain, acceptor_resi, name='O')
            cmd.distance(distance_name, donor, acceptor)
    cmd.deselect()

def assign_color(sse):
    """Assigns a color to, SSE based on 'color' field or generated from the label."""
    if COLOR in sse:
        return sse[COLOR]
    else:
        # magic_value = hash(sse[LABEL]) % 1000
        random.seed(sse[LABEL])
        magic_value = random.randint(0, 999)
        return 's' + str(magic_value).zfill(3)

def label_sse(sse, selection, center_of_mass=None):
    chain = sse[CHAIN_ID]
    start = sse[START_RESI]
    end = sse[END_RESI]
    middle = (start + end) // 2
    if center_of_mass is None:
        label_selection = '(' + selection + ') and chain ' + chain + ' and resi ' + str(middle) + ' and name CB'
        if cmd.count_atoms(label_selection)==0:
            label_selection = '(' + selection + ') and chain ' + chain + ' and resi ' + str(middle) + ' and name CA'
    else:
        CBs = '(' + selection + ') and chain ' + chain + ' and resi ' + resi_str(middle-1) + '-' + resi_str(middle+2) + ' and name CB'
        label_selection = get_farthest_atom(CBs, center_of_mass)
        if label_selection is None:
            CAs = '(' + selection + ') and chain ' + chain + ' and resi ' + resi_str(middle-1) + '-' + resi_str(middle+2) + ' and name CA'
            label_selection = get_farthest_atom(CAs, center_of_mass)
        if label_selection is None:
            label_selection = '(' + selection + ') and chain ' + chain + ' and resi ' + str(middle) + ' and name CA'
    cmd.label(label_selection, '"' + sse[LABEL].replace('"', '\\"') + '"')

def resi_str(number):
    """Convert residue number from int or str to PyMOL-style string, e.g. -5 -> '\-5'."""
    return str(number).replace('-', '\-')

def get_farthest_atom(candidate_selection, center_of_mass):
    x0, y0, z0 = center_of_mass
    namespace = {'coords': []}
    cmd.iterate_state(-1, candidate_selection, 'coords.append((resi, name, x, y, z))', space=namespace)
    sqdist_resi_name = [((x-x0)**2 + (y-y0)**2 + (z-z0)**2, resi, name) for (resi, name, x, y, z) in namespace['coords']]
    if len(sqdist_resi_name) > 0:
        sqdist, resi, name = max(sqdist_resi_name)
        return '(' + candidate_selection + ') and resi ' + resi + ' and name ' + name
    else:
        return None

def safe_object_name(name):
    """Replaces apostrophes, which are not allowed by PyMOL."""
    return name.replace('\'', SYMBOL_INSTEAD_APOSTROPHE)  # apostrophe is not allowed in PyMOL identifiers

def convert_ranges(ranges):
    """Converts changes from SecStrAnnotator format (1:10,15:20) to PyMOL format (1-10+15-20)."""
    return str(ranges).replace('-','\-').replace(':','-').replace(',','+')

def selection_expression(object_name, chain, ranges, symbol=None, name=None):
    """Creates PyMOL selection expression."""
    parts = ['('+object_name+')']
    if chain is not None:
        if chain=='.' or chain==' ':  # PyMOL ain't like dot and space values.
            chain = '""'
        parts.append('chain ' + chain)
    if ranges is not None:
        parts.append('resi ' + convert_ranges(ranges))
    if symbol is not None:
        parts.append('symbol ' + symbol)
    if name is not None:
        parts.append('name ' + name)
    return ' & '.join(parts) + ' '

def apply_representations(selection, protein_reprs=[], ligand_reprs=[], bonds_reprs=[]):
    """Hides everything of selection and shows given representations for HET and NOT HET atoms."""
    protein_sel = '(' + selection + ') & not het'
    ligand_sel = 'byres ((%s) expand %f & het & not resn HOH)' % (selection, LIGAND_EXTENSION)

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


def create_session(directory, template, query):
    """Creates and saves a session with superimposed template and query and selected and colored SSEs."""
    files = []

    if template is not None:
        template_id, template_chain, template_range = pdb_chain_ranges(template) 
        template_struct_file = path.join(directory, template_id + TEMPLATE_STRUCT_EXT)
        template_annot_file = path.join(directory, template_id + TEMPLATE_ANNOT_EXT)
        files.append(template_struct_file)
        files.append(template_annot_file)
    
    query_id, query_chain, query_range = pdb_chain_ranges(query) 
    query_struct_file = path.join(directory, query_id + QUERY_STRUCT_EXT)
    query_annot_file = path.join(directory, query_id + QUERY_ANNOT_EXT)
    session_file = path.join(directory, query_id + SESSION_EXT)
    template_png_file = path.join(directory, query_id + TEMPLATE_PNG_EXT)
    query_png_file = path.join(directory, query_id + QUERY_PNG_EXT)
    files.append(query_struct_file)
    files.append(query_annot_file)
    
    for filename in files:
        if not path.isfile(filename):
            print('ERROR: File not found: ' + filename)
            cmd.quit(1)
    
    query_object = QUERY_OBJECT_PREFIX + query_id
    query_annot_object = query_object + '.sses'

    if template is not None:
        template_object = TEMPLATE_OBJECT_PREFIX + template_id
        if template_object == query_object:
            template_object = template_object + '_'
        template_annot_object = template_object + '.sses'
        cmd.load(template_struct_file, template_object)
        if template_chain is not None:
            cmd.remove(f'{template_object} and not chain {template_chain}')
        template_selection = selection_expression(template_object, template_chain, None)
        color_by_annotation(template_id, template_selection, TEMPLATE_BASE_COLOR, template_annot_file, selection_prefix=template_annot_object)
        cmd.disable(template_object)
        cmd.disable(template_annot_object)

    cmd.load(query_struct_file, query_object)
    if query_chain is not None:
        cmd.remove(f'{query_object} and not chain {query_chain}')
    query_selection = selection_expression(query_object, query_chain, None)
    color_by_annotation(query_id, query_selection, QUERY_BASE_COLOR, query_annot_file, selection_prefix=query_annot_object)
    
    cmd.hide('everything', 'all')
    if template is not None:
        apply_representations(template_selection, protein_reprs=PROTEIN_REPRESENTATIONS, ligand_reprs=LIGAND_REPRESENTATIONS)
    apply_representations(query_selection, protein_reprs=PROTEIN_REPRESENTATIONS, ligand_reprs=LIGAND_REPRESENTATIONS, bonds_reprs=BOND_REPRESENTATIONS)
    cmd.dss()
    cmd.zoom('vis')
    cmd.save(session_file)
    if IMAGE_SIZE is not None:
        width, height = IMAGE_SIZE
        if template is not None:
            cmd.orient(template_selection)
        else:
            cmd.hide('labels')
            cmd.orient(query_selection)
        cmd.zoom('vis')
        cmd.ray(width, height)
        cmd.png(query_png_file)
        if template is not None:
            cmd.enable(template_object)
            cmd.enable(template_annot_object)
            cmd.disable(query_object)
            cmd.disable(query_annot_object)
            cmd.ray(width, height)
            cmd.png(template_png_file)


# Main script

try:
    if USE_CIF:
        cmd.set('cif_use_auth', False)
    create_session(directory, template, query)
    cmd.quit(0)
except Exception as e:
    print('ERROR: Exception raised: ' + str(e))
    cmd.quit(1)