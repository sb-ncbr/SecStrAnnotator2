import os
from os import path
import sys
import json
import requests
import argparse
from collections import defaultdict
from typing import Tuple

import no_gap_align

#  CONSTANTS  ##############################################################################

from constants import *

#  FUNCTIONS  ##############################################################################

class LazyDict:
    def __init__(self, initializer):
        self.initializer = initializer
        self.dictionary = {}
    def get(self, key):
        if key not in self.dictionary:
            self.dictionary[key] = self.initializer(key)
        return self.dictionary[key]
    def __getitem__(self, key):
        return self.get(key)

class Label2AuthConverter:
    def __init__(self, conversion_table_file):
        self.file = conversion_table_file
        self.table = {}
        with open(conversion_table_file) as r:
            for line in r:
                if not line.lstrip().startswith('#'):
                    chain, resi, auth_chain, auth_resi, auth_ins_code = line.strip().split('\t')
                    self.table[(chain, int(resi))] = (auth_chain, int(auth_resi), auth_ins_code)
        # sys.stderr.write(f'{conversion_table_file}:\n{self.table}\n')
    def auth_chain_resi_ins(self, chain, resi):
        # sys.stderr.write(f'{self.table}\n')
        try:
            return self.table[(chain, resi)]
        except KeyError:
            raise Exception(f'Did not find residue {chain} {resi} in {self.file}')

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('all_annotations_file', help='JSON file with the list of domains with annotations in SecStrAPI format')
parser.add_argument('reference_alignments_directory', help='Directory with reference alignments', type=str, default=None)
parser.add_argument('--labels', help='Comma-separated labels of SSEs to be processed (by default: all)', type=str, default=None)
parser.add_argument('--label2auth_dir', help='Directory with <PDB>.label2auth.tsv files for residue numbering conversion', type=str, default=None)
args = parser.parse_args()

all_annotations_file = args.all_annotations_file
alignments_dir = args.reference_alignments_directory
labels = args.labels.split(',') if args.labels is not None else None
label2auth_dir = args.label2auth_dir

#  MAIN  ##############################################################################

with open(all_annotations_file) as r:
    all_annotations = json.load(r)

aligner = no_gap_align.NoGapAligner()
realigners_manager = LazyDict(lambda label: no_gap_align.Realigner(path.join(alignments_dir, label + '.fasta')))

api_version = all_annotations[API_VERSION]
pdb2domains = all_annotations[ANNOTATIONS]
for pdb, domains in pdb2domains.items():
    dom_list = domains.values() if isinstance(domains, dict) else domains[:]
    converter = Label2AuthConverter(path.join(label2auth_dir, pdb + '.label2auth.tsv')) if label2auth_dir is not None else None
    for domain in dom_list:
        for sse in domain[SSES]:
            label = sse.get(LABEL, None)
            if label is not None and (labels is None or label in labels):
                realigner = realigners_manager[label]
                shift, pivot_index = realigner.aligning_shift_and_pivot(sse[SEQUENCE])
                pivot_residue = sse[START] + pivot_index
                if not sse[START] <= pivot_residue <= sse[END]:
                    message = f'{pdb}: pivot residue of {label} ({pivot_residue}) falls out out {label} ({sse[START]}-{sse[END]})'
                    # raise Exception(message)
                    sys.stderr.write(f'WARNING: {message}\n')
                    # continue
                sse[PIVOT_RESIDUE] = pivot_residue
                if converter is not None:
                    _, sse[AUTH_PIVOT_RESIDUE], sse[AUTH_PIVOT_RESIDUE_INS_CODE] = converter.auth_chain_resi_ins(sse[SSE_CHAIN], sse[PIVOT_RESIDUE])

json.dump(all_annotations, sys.stdout, indent=4)
print()

n_pdbs = len(all_annotations[ANNOTATIONS])
n_domains = sum( len(doms) for doms in all_annotations[ANNOTATIONS].values() )
n_labels = len(realigners_manager.dictionary)
sys.stderr.write(f'Added pivot residues for {n_domains} domains in {n_pdbs} PDB entries ({n_labels} labels)\n')