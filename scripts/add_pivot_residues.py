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

class RealignersManager:
    def __init__(self, alignments_directory):
        self.alignments_directory = alignments_directory
        self.label2realigner = {}
    def get_realigner(self, label) -> no_gap_align.Realigner:
        if label not in self.label2realigner:
            self.label2realigner[label] = no_gap_align.Realigner(path.join(self.alignments_directory, label + '.fasta'))
            pivot = self.label2realigner[label].pivot_index
        return self.label2realigner[label]

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('all_annotations_file', help='JSON file with the list of domains with annotations in SecStrAPI format')
parser.add_argument('reference_alignments_directory', help='Directory with reference alignments', type=str, default=None)
parser.add_argument('--labels', help='Comma-separated labels of SSEs to be processed (by default: all)', type=str, default=None)
args = parser.parse_args()

all_annotations_file = args.all_annotations_file
alignments_dir = args.reference_alignments_directory
labels = args.labels.split(',') if args.labels is not None else None

#  MAIN  ##############################################################################

with open(all_annotations_file) as r:
    all_annotations = json.load(r)

aligner = no_gap_align.NoGapAligner()
realigners_manager = RealignersManager(alignments_dir)

api_version = all_annotations[API_VERSION]
pdb2domains = all_annotations[ANNOTATIONS]
for pdb, domains in pdb2domains.items():
    dom_list = domains.values() if isinstance(domains, dict) else domains[:]
    for domain in dom_list:
        for sse in domain[SSES]:
            label = sse.get(LABEL, None)
            if label is not None and (labels is None or label in labels):
                realigner = realigners_manager.get_realigner(label)
                shift, pivot_index = realigner.aligning_shift_and_pivot(sse[SEQUENCE])
                sse[PIVOT_RESIDUE] = sse[START] + pivot_index
                # if label == 'A':
                #     sys.stderr.write(f'{pdb} {label} {shift} {pivot_index} ({sse[START]} - {sse[PIVOT_RESIDUE]} - {sse[END]})\n')

json.dump(all_annotations, sys.stdout, indent=4)
print()

n_pdbs = len(all_annotations[ANNOTATIONS])
n_domains = sum( len(doms) for doms in all_annotations[ANNOTATIONS].values() )
n_labels = len(realigners_manager.label2realigner)
sys.stderr.write(f'Added pivot residues for {n_domains} domains in {n_pdbs} PDB entries ({n_labels} labels)\n')