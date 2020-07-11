'''
This Python3 script aligns SSE sequences to reference alignments and adds generic numbering information into the annotation file.

Example usage:
    python3  add_reference_residues.py  annotations.json  aligments/  --labels A,B,C  --label2auth_dir structures/
'''

import argparse
from typing import Dict, Any, Optional, Tuple, List, Union
import os
from os import path
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
    parser.add_argument('reference_alignments_dir', help='Directory with reference alignments', type=str, default=None)
    parser.add_argument('--labels', help='Comma-separated labels of SSEs to be processed (by default: all)', type=str, default=None)
    parser.add_argument('--label2auth_dir', help='Directory with <PDB>.label2auth.tsv files for residue numbering conversion', type=str, default=None)
    args = parser.parse_args()
    return vars(args)


def main(all_annotations_file: str, reference_alignments_dir: str, labels: Union[str, List[str], None] = None, label2auth_dir: str = None) -> Optional[int]:
    '''Align SSE sequences to reference alignments and add generic numbering information into the annotation file.'''

    with open(all_annotations_file, 'r', encoding=lib.DEFAULT_ENCODING) as r:
        all_annotations = json.load(r)

    do_all_labels = labels is None or labels == 'all'
    if not do_all_labels and isinstance(labels, str):
        labels = labels.split(',')

    aligner = no_gap_align.NoGapAligner()
    realigners_manager = lib.LazyDict(lambda label: no_gap_align.Realigner(path.join(reference_alignments_dir, label + '.fasta')))

    api_version = all_annotations[API_VERSION]
    pdb2domains = all_annotations[ANNOTATIONS]
    for pdb, domains in pdb2domains.items():
        dom_list = domains.values() if isinstance(domains, dict) else domains[:]
        converter = lib.Label2AuthConverter(path.join(label2auth_dir, pdb + '.label2auth.tsv'), unknown_ins_code_as_empty_string=True) if label2auth_dir is not None else None
        for domain in dom_list:
            for sse in domain[SSES]:
                label = sse.get(LABEL, None)
                if label is not None and (do_all_labels or label in labels):
                    realigner = realigners_manager[label]
                    shift, ref_index = realigner.aligning_shift_and_pivot(sse[SEQUENCE])
                    ref_residue = sse[START] + ref_index
                    if not sse[START] <= ref_residue <= sse[END]:
                        message = f'{pdb}: reference residue of {label} ({ref_residue}) falls out of {label} ({sse[START]}-{sse[END]})'
                        sys.stderr.write(f'    WARNING: {message}\n')
                    sse[PIVOT_RESIDUE] = ref_residue
                    if converter is not None:
                        try:
                            _, sse[AUTH_PIVOT_RESIDUE], sse[AUTH_PIVOT_RESIDUE_INS_CODE] = converter.auth_chain_resi_ins(sse[CHAIN_ID], sse[PIVOT_RESIDUE])
                        except KeyError:
                            message = f'{pdb}: reference residue of {label} ({ref_residue}) is not modelled in the structure)'
                            sys.stderr.write(f'  WARNING: {message}\n')

    json.dump(all_annotations, sys.stdout, indent=4)
    print()

    n_pdbs = len(all_annotations[ANNOTATIONS])
    n_domains = sum( len(doms) for doms in all_annotations[ANNOTATIONS].values() )
    n_labels = len(realigners_manager.dictionary)
    sys.stderr.write(f'Added reference residues for {n_domains} domains in {n_pdbs} PDB entries ({n_labels} labels)\n')


if __name__ == '__main__':
    args = parse_args()
    exit_code = main(**args)
    if exit_code is not None:
        exit(exit_code)