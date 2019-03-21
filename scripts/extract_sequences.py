import os
import sys
import json
import requests
import argparse
from collections import defaultdict
from typing import Tuple
from os import path
import shutil

#  CONSTANTS  ##############################################################################

from constants import *

SEQUENCE_EXT = '.fasta'

#  FUNCTIONS  ##############################################################################

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('all_annotations_file', help='JSON file with the list of domains with annotations in SecStrAPI format', type=str)
parser.add_argument('output_directory', help='Directory to output sequences (a multi-FASTA per label)', type=str)
args = parser.parse_args()

all_annotations_file = args.all_annotations_file
output_directory = args.output_directory

#  MAIN  ##############################################################################

with open(all_annotations_file) as r:
    all_annotations = json.load(r)

shutil.rmtree(output_directory, ignore_errors=True)
os.makedirs(output_directory)

label2domain_sequence = defaultdict(lambda: [])

for pdb, domains in all_annotations[ANNOTATIONS].items():
    dom_list = domains.values() if isinstance(domains, dict) else domains[:]
    for domain in dom_list:
        name = ','.join((domain[PDB], domain[CHAIN], domain[RANGES]))
        for sse in domain[SSES]:
            if LABEL in sse and sse[LABEL] is not None:
                label2domain_sequence[sse[LABEL]].append((name, sse[SEQUENCE]))

for label, domnames_sequences in label2domain_sequence.items():
    with open(path.join(output_directory, label + SEQUENCE_EXT), 'w') as w:
        for domname, sequence in domnames_sequences:
            w.write(f'>{domname}\n{sequence}\n')

n_pdbs = len(all_annotations[ANNOTATIONS])
n_domains = sum( len(doms) for doms in all_annotations[ANNOTATIONS].values() )
n_labels = len(label2domain_sequence)
sys.stderr.write(f'Extracted sequences for {n_domains} domains in {n_pdbs} PDB entries ({n_labels} labels)\n')