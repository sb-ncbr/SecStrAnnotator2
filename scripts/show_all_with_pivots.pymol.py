import json
import os
from os import path
import sys
from pymol import cmd

import annotate_sec_str_extension as secstr

cmd.set('cif_use_auth', False)

# How to run: pymol  -qcr show_all_with_pivots.pymol.py  --  ANNOTATION_FILE  STRUCTURE_DIR  CONSENSUS  OUT_SESSION

if len(sys.argv) != 5:
	print('Usage: pymol  -qcr show_all_with_pivots.pymol.py  --  ANNOTATION_FILE  STRUCTURE_DIR  CONSENSUS  OUT_SESSION')
	print(len(sys.argv))
	exit()
_, annotation_file, structure_dir, consensus_file, out_session = sys.argv

with open(annotation_file) as f:
	js = json.load(f)
	annotations = js['annotations']

pdbs = list(annotations.keys())

# cmd.load(consensus_file, 'consensus')

for pdb, domains in sorted(annotations.items()):
	# print(pdb)
	# continue
	if len(domains) != 1:
		raise NotImplementedError(pdb)
	domain = domains[0]
	domdef = pdb + ',' + domain[secstr.CHAIN] + ',:-aligned'
	structure_file = path.join(structure_dir, domdef + '.cif')
	structure = pdb
	domain_name = domain[secstr.PDB] + domain[secstr.CHAIN]  #TODO read domain_name
	print(domain_name)
	cmd.load(structure_file, structure)
	secstr.annotate_sec_str(structure, annotation_file=annotation_file, name=pdb)  #TODO read per-pdb annotation file
	# cmd.cealign('consensus', domain_name)
cmd.save(out_session)