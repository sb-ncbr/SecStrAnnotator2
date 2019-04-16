import json
import os
from os import path
import sys
from pymol import cmd
cmd.set('cif_use_auth', False)

# How to run: pymol  -qcr show_all_with_pivots.pymol.py  --  ANNOTATION_FILE  STRUCTURE_DIR  CONSENSUS  OUT_SESSION

if len(sys.argv) != 5:
	print('Usage: pymol  -qcr show_all_with_pivots.pymol.py  --  ANNOTATION_FILE  STRUCTURE_DIR  CONSENSUS  OUT_SESSION')
	print(len(sys.argv))
	exit()
_, annotation_file, structure_dir, consensus_file, out_session = sys.argv

with open(annotation_file) as f:
	annotations = json.load(f)

pdbs = list(annotations['annotations'].keys())

cmd.load(consensus_file, 'consensus')

for pdb, domain, chain, rang in domains:
	print(domain)  # debug
	in_file = path.join(in_directory, domain + '.cif')
	out_file = path.join(out_directory, domain + '.cif')
	cmd.load(in_file, 'protein')
	cmd.cealign('consensus', 'protein')
	cmd.save(out_file, 'protein')
	cmd.delete('protein')

cmd.quit()