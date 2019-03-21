import os
from os import path
import shutil
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

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('all_annotations_file', help='JSON file with the list of domains with annotations in SecStrAPI format')
parser.add_argument('--labels', help='Comma-separated labels of SSEs to be processed (by default: all)', type=str, default=None)
parser.add_argument('--alignments_directory', help='Directory to output alignments', type=str, default=None)
parser.add_argument('--trees_directory', help='Directory to output trees', type=str, default=None)
parser.add_argument('--logos_directory', help='Directory to output sequence logos', type=str, default=None)
args = parser.parse_args()

all_annotations_file = args.all_annotations_file
labels = args.labels.split(',') if args.labels is not None else None
alignments_dir = args.alignments_directory
trees_dir = args.trees_directory
logos_dir = args.logos_directory

#  MAIN  ##############################################################################

with open(all_annotations_file) as r:
    all_annotations = json.load(r)

pdb2domains = all_annotations[ANNOTATIONS]
domains = [ dom for doms in pdb2domains.values() for dom in doms ]
label2seqs = defaultdict(lambda: [])
for domain in domains:
    for sse in domain[SSES]:
        name = domain[PDB] + ',' + domain[CHAIN] + ',' + domain[RANGES]
        label2seqs[sse[LABEL]].append((name, sse[SEQUENCE]))

if labels is None:
    labels = sorted(label2seqs.keys())

if alignments_dir is not None:
    shutil.rmtree(alignments_dir, ignore_errors=True)
    os.makedirs(alignments_dir)
if trees_dir is not None:
    shutil.rmtree(trees_dir, ignore_errors=True)
    os.makedirs(trees_dir)
if logos_dir is not None:
    shutil.rmtree(logos_dir, ignore_errors=True)
    os.makedirs(logos_dir)

aligner = no_gap_align.NoGapAligner()
for label in labels:
    print(label)
    names, sequences = zip(*label2seqs[label])
    # for name, sequence in zip(names, sequences):
    #     if '?' in sequence:
    #         print(name, sequence)
    aligner.align(sequences, names=names)
    # aligner.print_tree()
    if alignments_dir is not None:
        aligner.output_alignment(path.join(alignments_dir, label + '.fasta'))
    if alignments_dir is not None:
        aligner.print_tree(output_file=path.join(trees_dir, label + '.txt'))
    if logos_dir is not None:
        aligner.output_logo(path.join(logos_dir, label + '.png'))
