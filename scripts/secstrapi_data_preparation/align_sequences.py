'''
This Python3 script reads SSE annotations in SecStrAPI format and performs multiple sequence alignment of SSE sequences (separately for each SSE label).

Example usage:
    python3  align_sequences.py annotations.json  --labels A,B,C  --alignments_dir aligments/  --trees_dir trees_NR/  --logos_dir logos_NR/
'''

import argparse
from typing import Dict, Any, Optional, Tuple, List, Union
import os
from os import path
import shutil
import sys
import json
from collections import defaultdict

import lib
import no_gap_align
from constants import *

#  CONSTANTS  ################################################################################


#  FUNCTIONS  ################################################################################


#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('all_annotations_file', help='JSON file with the list of domains with annotations in SecStrAPI format')
    parser.add_argument('--labels', help='Comma-separated labels of SSEs to be processed (default: all)', type=str, default=None)
    parser.add_argument('--alignments_dir', help='Directory to output alignments', type=str, default=None)
    parser.add_argument('--trees_dir', help='Directory to output trees', type=str, default=None)
    parser.add_argument('--logos_dir', help='Directory to output sequence logos', type=str, default=None)
    parser.add_argument('--matrices_dir', help='Directory to output alignment matrices', type=str, default=None)
    parser.add_argument('--labels_for_matrices', help='Comma-separated labels of SSEs to create matrix files for (default: =labels)', type=str, default=None)
    parser.add_argument('--ref_residue', help='Number assigned to the reference residue (i.e. the most conserved), default: 50', type=int, default=50)
    args = parser.parse_args()
    return vars(args)


def main(all_annotations_file: str, labels: Union[str, List[str], None] = None, alignments_dir: Optional[str] = None, 
        trees_dir: Optional[str] = None, logos_dir: Optional[str] = None, matrices_dir: Optional[str] = None, labels_for_matrices: Union[str, List[str], None] = None, ref_residue: int = 50) -> Optional[int]:
    '''Read SSE annotations in SecStrAPI format and perform multiple sequence alignment of SSE sequences (separately for each SSE label).'''

    LOGO_UNITS = 'bits'

    with open(all_annotations_file, 'r', encoding=lib.DEFAULT_ENCODING) as r:
        all_annotations = json.load(r)

    pdb2domains = all_annotations[ANNOTATIONS]
    label2seqs = defaultdict(list)
    for pdb, domains in pdb2domains.items():
        for name, domain in lib.iterate_names_domains(domains):
            for sse in domain[SSES]:
                label2seqs[sse[LABEL]].append((name, sse[SEQUENCE]))

    if labels is None:
        labels = sorted(label2seqs.keys())
    elif isinstance(labels, str):
        labels = labels.split(',')

    if labels_for_matrices is None:
        labels_for_matrices = labels
    elif isinstance(labels_for_matrices, str):
        labels_for_matrices = labels_for_matrices.split(',')

    if alignments_dir is not None:
        shutil.rmtree(alignments_dir, ignore_errors=True)
        os.makedirs(alignments_dir)
    if trees_dir is not None:
        shutil.rmtree(trees_dir, ignore_errors=True)
        os.makedirs(trees_dir)
    if logos_dir is not None:
        shutil.rmtree(logos_dir, ignore_errors=True)
        os.makedirs(logos_dir)
    if matrices_dir is not None:
        shutil.rmtree(matrices_dir, ignore_errors=True)
        os.makedirs(matrices_dir)

    aligner = no_gap_align.NoGapAligner()
    aligner.print_column_statistics(only_header=True)
    for label in labels:
        # print(label, file=sys.stderr)
        names, sequences = zip(*label2seqs[label])
        aligner.align(sequences, names=names)
        aligner.print_column_statistics(label=label)
        if alignments_dir is not None:
            aligner.output_alignment(path.join(alignments_dir, label + '.fasta'))
        if trees_dir is not None:
            aligner.print_tree(output_file=path.join(trees_dir, label + '.txt'))
        if logos_dir is not None:
            aligner.output_logo(path.join(logos_dir, label + '.png'), tool='logomaker', pivot_as=ref_residue, units=LOGO_UNITS)
            aligner.output_logo(path.join(logos_dir, label + '.tif'), tool='logomaker', pivot_as=ref_residue, units=LOGO_UNITS)
        if matrices_dir is not None:
            aligner.output_alignment_matrices(path.join(matrices_dir, label + '.json'), pivot_as=ref_residue, units=LOGO_UNITS)
    if matrices_dir is not None:
        all_matrices = {}
        for label in labels_for_matrices:
            with open(path.join(matrices_dir, label + '.json')) as r:
                matrices = json.load(r)
            all_matrices[label] = matrices
        with open(path.join(matrices_dir, 'ALL.json'), 'w') as w:
            json.dump(all_matrices, w, separators=(',',':'))


if __name__ == '__main__':
    args = parse_args()
    exit_code = main(**args)
    if exit_code is not None:
        exit(exit_code)