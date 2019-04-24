import json
import os
from os import path
import sys
from pymol import cmd

import annotate_sec_str_extension as secstr

cmd.set('cif_use_auth', False)

# How to run: pymol  -qcr show_all_with_pivots.pymol.py  --  ANNOTATION_FILE  STRUCTURE_DIR  CONSENSUS  OUT_SESSION
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

with open(annotation_file) as f:
	js = json.load(f)
	annotations = js['annotations']

# pdbs = list(annotations.keys())
n = len(annotations)

# cmd.load(consensus_file, 'consensus')

for i, (pdb, domains) in enumerate(sorted(annotations.items())):
	if len(domains) != 1:
		raise NotImplementedError(pdb)
	domain_name, domain = domains.items()[0]
	structure_file = path.join(structure_dir, domain_name + '-aligned.cif')
	structure = pdb
	print(str(100*i//n) + '%  ' + domain_name)
	cmd.load(structure_file, structure)
	secstr.annotate_sec_str(structure, annotation_file=path.join(separate_annotations_dir, pdb + '.json'), name=domain_name)
	# cmd.cealign('consensus', domain_name)
	cmd.hide('lines', structure)
	cmd.hide('nonbonded', structure)
	cmd.show('ribbon', domain_name)

cmd.set('ribbon_width', 0.5)
cmd.save(out_session)