# PymolScriptAlign for SecStrAnnotator 1.0
# Superimposes template and query domain.
# Usage: 
#     pymol -qcyr script_align.py -- file_format align_method directory template_pdbid template_chain template_ranges query_pdbid query_chain query_ranges


# Imports

from os import path 
import sys
from pymol import cmd


# Get command line arguments

N_PSEUDOARGUMENTS = 1
arguments = sys.argv[N_PSEUDOARGUMENTS:]
if len(arguments) != 9:
	print('ERROR: Exactly 9 command line arguments required, ' + str(len(arguments)) + ' given: ' + ' '.join(arguments))
	cmd.quit(1)
file_format, method, directory, template_id, template_chain, template_range, query_id, query_chain, query_range = arguments


# Additional constants

if file_format == 'cif':
	USE_CIF = True
elif file_format == 'pdb':
	USE_CIF = False
else:
	raise Exception('Unknown file format: ' + file_format)

STRUCTURE_EXT = '.cif' if USE_CIF else '.pdb'
ALIGNED_STRUCTURE_EXT = '-aligned.cif' if USE_CIF else '-aligned.pdb'
ALIGNMENT_OBJECT = 'aln'


# Auxiliary functions

def convert_ranges(ranges):
	"""Converts changes from SecStrAnnotator format (1:10,15:20) to PyMOL format (1-10+15-20)."""
	return ranges.replace('-','\-').replace(':','-').replace(',','+')

def selection_expression(object_name, chain, ranges, symbol=None):
	"""Creates PyMOL selection expression."""
	parts = [object_name, 'chain ' + chain, 'resi ' + convert_ranges(ranges)]
	if symbol is not None:
		parts.append('symbol ' + symbol)
	# print(' & '.join(parts) + ' ')
	return ' & '.join(parts) + ' '

def aligned_residues_from_alignment_object(alignment_object_name):
	"""Returns a list of paired residues as 4-tuples (template_chain, template_resi, query_chain, query_resi)."""
	raw_aln = cmd.get_raw_alignment(alignment_object_name)
	cr_aln = []
	for idx1, idx2 in raw_aln:
		cmd.iterate(idx1[0]+' & index '+str(idx1[1]), PASSED_VAR+'=(chain,resi)')
		crt = stored.passed_var
		cmd.iterate(idx2[0]+' & index '+str(idx2[1]), PASSED_VAR+'=(chain,resi)')
		crq = stored.passed_var
		cr_aln.append((crt[0], int(crt[1]), crq[0], int(crq[1])))
	return cr_aln

def perform_alignment(method, directory, template_id, template_chain, template_range, query_id, query_chain, query_range):
	"""Perform specified structural alignment method."""
	t = 'template'  # name for whole template
	q = 'query'  # name for whole query protein
	template_struct_file = path.join(directory, template_id + STRUCTURE_EXT)
	query_struct_file = path.join(directory, query_id + STRUCTURE_EXT)
	query_aligned_struct_file = path.join(directory, query_id + ALIGNED_STRUCTURE_EXT)
	for filename in [template_struct_file, query_struct_file]:
		if not path.isfile(filename):
			print('ERROR: File not found: ' + filename)
			cmd.quit(1)
	cmd.load(template_struct_file, t)
	cmd.load(query_struct_file, q)

	sel_t = selection_expression(t, template_chain, template_range)
	sel_q = selection_expression(q, query_chain, query_range)

	if cmd.count_atoms(sel_t) == 0:
		print('ERROR: Template domain contains no atoms.')
		cmd.quit(1)
	if cmd.count_atoms(sel_q) == 0:
		print('ERROR: Query domain contains no atoms.')
		cmd.quit(1)

	if method == 'align':
		cmd.align(sel_q, sel_t)
	elif method == 'super':
		cmd.super(sel_q, sel_t)
	elif method == 'cealign':
		aln_result = cmd.cealign(sel_t, sel_q, object=ALIGNMENT_OBJECT)
	else:
		print('ERROR: Unknown structural alignment method: '+method)
		cmd.quit(1)

	cmd.save(query_aligned_struct_file, q)


# Main script

try:
	if USE_CIF:
		cmd.set('cif_use_auth', False)
	perform_alignment(method, directory, template_id, template_chain, template_range, query_id, query_chain, query_range)
	cmd.quit(0)
except Exception as e:
	print('ERROR: Exception raised: ' + str(e))
	cmd.quit(1)