import json
import os
from os import path
import sys
from pymol import cmd

import annotate_sec_str_extension as secstr

cmd.set('cif_use_auth', False)

# How to run: pymol  -qcr show_all_with_pivots.pymol.py  --  DIRECTORY  OUT_SESSION
argument_names = 'DIRECTORY  OUT_SESSION'

if len(sys.argv) != 1 + len(argument_names.split()):
	print('Usage: pymol  -qcr show_all_with_pivots.pymol.py  --  ' + argument_names)
	print(len(sys.argv))
	exit()
# _, annotation_file, structure_dir, consensus_file, out_session = sys.argv
_, directory, out_session = sys.argv

annotation_file = path.join(directory, 'annotations_with_pivots_best.json')
structure_dir = path.join(directory, 'structures')
separate_annotations_dir = path.join(directory, 'annotations_all')
consensus_file = path.join(directory, 'consensus.cif')

def single(iterable):
	iterator = iter(iterable)
	try:
		value = next(iterator)
		try:
			_ = next(iterator)
			raise ValueError('Iterable contains more than one element)')
		except StopIteration:
			return value
	except StopIteration:
		raise ValueError('Iterable contains no elements')


with open(annotation_file) as f:
	js = json.load(f)
	annotations = js['annotations']

# pdbs = list(annotations.keys())
n = len(annotations)

SEPARATE = False

if SEPARATE:
	extra_residues_before = 25
	extra_residues_after = 25
	qualities = {}
	with open(path.join(directory, 'qualities.tsv')) as r:
		for line in r:
			pdb, quality = line.strip().split()
			qualities[pdb] = float(quality)
	labels = { sse['label'] for domains in annotations.values() for domain in domains.values() for sse in domain['secondary_structure_elements'] }
	labels = "A B''' C G H".split()  # for after-align
	labels = "B B'' E F".split()  # for before-align

	for i, (pdb, domains) in enumerate(sorted(annotations.items())):
		domain_name, domain = single(domains.items())
		structure_file = path.join(structure_dir, domain_name + '-aligned.cif')
		structure = pdb
		print(str(100*i//n) + '%  ' + domain_name)
		cmd.load(structure_file, structure)
		domain_selection = pdb + ' and chain ' + domain['chain_id']
		cmd.select(domain_name, domain_selection)
		cmd.hide('everything', structure)
		cmd.show('ribbon', domain_name)
	
	for the_label in sorted(labels):
		print(the_label)
		occurring = [ (domain_name, domain) for domains in annotations.values() for domain_name, domain in domains.items() if any( sse['label']==the_label for sse in domain['secondary_structure_elements'] ) ]
		best_name, best_domain = max(occurring, key=lambda nd: qualities[nd[0][0:4]])
		print(best_name)
		cmd.hide('sticks')
		cmd.delete('*.*')
		sse = single( sse for sse in best_domain['secondary_structure_elements'] if sse['label']==the_label  )
		ref_selection = best_name + ' and resi ' + str(sse['start']-extra_residues_before) + '-' + str(sse['end']+extra_residues_after)
		cmd.select('ref', ref_selection)
		for i, (pdb, domains) in enumerate(sorted(annotations.items())):
			domain_name, domain = single(domains.items())
			if any( sse['label']==the_label for sse in domain['secondary_structure_elements'] ):
				print(str(100*i//n) + '%  ' + domain_name)
				domain_selection = pdb + ' and chain ' + domain['chain_id']
				sse = single( sse for sse in domain['secondary_structure_elements'] if sse['label']==the_label )
				sse_name = domain_name + '.' + the_label.replace("'", "+")
				sse_selection = domain_name + ' and resi ' + str(sse['start']) + '-' + str(sse['end'])
				cmd.select(sse_name, sse_selection)
				cmd.color('gray80', domain_name)
				cmd.color('green', sse_name + ' and symbol C')
				if 'pivot_residue' in sse:
					pivot = sse['pivot_residue']
					pivot_name = domain_name + '.' + the_label.replace("'", "+") + '.pivot'
					pivot_selection = domain_name + ' and resi ' + str(sse['pivot_residue'])
					cmd.select(pivot_name, pivot_selection)
					cmd.show('sticks', pivot_name + ' and  name CA+CB')
					cmd.color('red', pivot_name + ' and name CB')
				mobile_selection = domain_name + ' and resi ' + str(sse['start']-extra_residues_before) + '-' + str(sse['end']+extra_residues_after)
				cmd.cealign('ref', mobile_selection)
		cmd.zoom('*.' + the_label.replace("'", "+"))
		cmd.deselect()
		cmd.set('ribbon_width', 0.5)
		cmd.save(out_session + '-' + the_label + '.pse')


# if not SEPARATE:
# 	labels = { sse['label'] for domains in annotations.values() for domain in domains.values() for sse in domain['secondary_structure_elements'] }
# 	# labels = "A B''' C G H".split()  # for after-align
# 	# labels = "B B'' E F".split()  # for before-align

# 	cmd.load(consensus_file, 'consensus')

# 	for i, (pdb, domains) in enumerate(sorted(annotations.items())):
# 		domain_name, domain = single(domains.items())
# 		structure_file = path.join(structure_dir, domain_name + '-aligned.cif')
# 		structure = pdb
# 		print(str(100*i//n) + '%  ' + domain_name)
# 		cmd.load(structure_file, structure)
# 		domain_selection = pdb + ' and chain ' + domain['chain_id']
# 		cmd.select(domain_name, domain_selection)
# 		cmd.cealign('consensus', domain_name)
# 		cmd.hide('everything', structure)
# 		cmd.show('ribbon', domain_name)
	
# 	for the_label in sorted(labels):
# 		print(the_label)
# 		occurring = [ (domain_name, domain) for domains in annotations.values() for domain_name, domain in domains.items() if any( sse['label']==the_label for sse in domain['secondary_structure_elements'] ) ]
# 		best_name, best_domain = max(occurring, key=lambda nd: qualities[nd[0][0:4]])
# 		print(best_name)
# 		cmd.hide('sticks')
# 		cmd.delete('*.*')
# 		sse = single( sse for sse in best_domain['secondary_structure_elements'] if sse['label']==the_label  )
# 		ref_selection = best_name + ' and resi ' + str(sse['start']-extra_residues_before) + '-' + str(sse['end']+extra_residues_after)
# 		cmd.select('ref', ref_selection)
# 		for i, (pdb, domains) in enumerate(sorted(annotations.items())):
# 			domain_name, domain = single(domains.items())
# 			if any( sse['label']==the_label for sse in domain['secondary_structure_elements'] ):
# 				print(str(100*i//n) + '%  ' + domain_name)
# 				domain_selection = pdb + ' and chain ' + domain['chain_id']
# 				sse = single( sse for sse in domain['secondary_structure_elements'] if sse['label']==the_label )
# 				sse_name = domain_name + '.' + the_label.replace("'", "+")
# 				sse_selection = domain_name + ' and resi ' + str(sse['start']) + '-' + str(sse['end'])
# 				cmd.select(sse_name, sse_selection)
# 				cmd.color('gray80', domain_name)
# 				cmd.color('green', sse_name + ' and symbol C')
# 				if 'pivot_residue' in sse:
# 					pivot = sse['pivot_residue']
# 					pivot_name = domain_name + '.' + the_label.replace("'", "+") + '.pivot'
# 					pivot_selection = domain_name + ' and resi ' + str(sse['pivot_residue'])
# 					cmd.select(pivot_name, pivot_selection)
# 					cmd.show('sticks', pivot_name + ' and  name CA+CB')
# 					cmd.color('red', pivot_name + ' and name CB')
# 				mobile_selection = domain_name + ' and resi ' + str(sse['start']-extra_residues_before) + '-' + str(sse['end']+extra_residues_after)
# 				cmd.cealign('ref', mobile_selection)
# 		cmd.zoom('*.' + the_label.replace("'", "+"))
# 		cmd.deselect()
# 		cmd.set('ribbon_width', 0.5)
# 		cmd.save(out_session + '-' + the_label + '.pse')



if not SEPARATE:
	cmd.load(consensus_file, 'consensus')
	for i, (pdb, domains) in enumerate(sorted(annotations.items())):
		if len(domains) != 1:
			raise NotImplementedError(pdb)
		domain_name, domain = domains.items()[0]
		structure_file = path.join(structure_dir, domain_name + '-aligned.cif')
		structure = pdb
		print(str(100*i//n) + '%  ' + domain_name)
		cmd.load(structure_file, structure)
		secstr.annotate_sec_str(structure, annotation_file=path.join(separate_annotations_dir, pdb + '.json'), name=domain_name)
		cmd.cealign('consensus', domain_name)
		cmd.hide('lines', structure)
		cmd.hide('nonbonded', structure)
		cmd.show('ribbon', domain_name)
	cmd.set('ribbon_width', 0.5)
	cmd.set('sphere_scale', 0.2)
	cmd.save(out_session)