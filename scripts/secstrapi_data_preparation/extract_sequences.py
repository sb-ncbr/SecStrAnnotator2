'''
This Python3 script extracts the amino acid sequences of annotated SSEs and prints them in FASTA format (one file per SSE label).

Example usage:
    python3  extract_sequences.py  annotations_ALL.json  sequences_ALL/
'''

import argparse
from typing import Dict, Any, Optional
import os
import sys
import json
from collections import defaultdict
from os import path
import shutil

import lib
from constants import *

#  CONSTANTS  ################################################################################

SEQUENCE_EXT = '.fasta'

#  FUNCTIONS  ################################################################################


#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('all_annotations_file', help='JSON file with the list of domains with annotations in SecStrAPI format', type=str)
    parser.add_argument('output_directory', help='Directory to output sequences (a multi-FASTA per label)', type=str)
    args = parser.parse_args()
    return vars(args)


def main(all_annotations_file: str, output_directory: str) -> Optional[int]:
    '''Extract the amino acid sequences of annotated SSEs and print them in FASTA format (one file per SSE label).'''

    with open(all_annotations_file, 'r', encoding=lib.DEFAULT_ENCODING) as r:
        all_annotations = json.load(r)

    shutil.rmtree(output_directory, ignore_errors=True)
    os.makedirs(output_directory)

    label2domain_sequence = defaultdict(lambda: [])

    for pdb, domains in all_annotations[ANNOTATIONS].items():
        for name, domain in lib.iterate_names_domains(domains):
            for sse in domain[SSES]:
                if LABEL in sse and sse[LABEL] is not None:
                    label2domain_sequence[sse[LABEL]].append((name, sse[SEQUENCE]))

    for label, domnames_sequences in label2domain_sequence.items():
        with open(path.join(output_directory, label + SEQUENCE_EXT), 'w', encoding=lib.DEFAULT_ENCODING) as w:
            for domname, sequence in domnames_sequences:
                w.write(f'>{domname}\n{sequence}\n')

    n_pdbs = len(all_annotations[ANNOTATIONS])
    n_domains = sum( len(doms) for doms in all_annotations[ANNOTATIONS].values() )
    n_labels = len(label2domain_sequence)
    sys.stderr.write(f'Extracted sequences for {n_domains} domains in {n_pdbs} PDB entries ({n_labels} labels)\n')


if __name__ == '__main__':
    args = parse_args()
    exit_code = main(**args)
    if exit_code is not None:
        exit(exit_code)