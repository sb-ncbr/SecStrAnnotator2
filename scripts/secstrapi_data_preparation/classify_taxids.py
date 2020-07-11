'''
This Python3 script reads a table domain|taxid and adds column with superkingdom (Arch/Bact/Euka/Viru), according to NCBI taxonomy hierarchy file (nodes.dmp).

Example usage:
    python3  classify_taxids.py  nodes.dmp  domains_taxids.tsv
'''

import argparse
from typing import Dict, Any
import os
import sys
import json

import lib

#  CONSTANTS  ################################################################################

EUKARYOTA = '2759'
BACTERIA = '2'
ARCHEA = '2157'
HUMAN = '9606'
VIRUSES = '10239'
GROUPS = [(EUKARYOTA,'Euka'), (BACTERIA,'Bact'), (ARCHEA,'Arch'), (VIRUSES, 'Viru')]

#  FUNCTIONS  ################################################################################

def parse_taxonomy_file(f):
    parent_dict = {}
    for line in iter(f.readline, ''):
        node, parent, *_ = (x.strip() for x in line.split('|'))
        parent_dict[node] = parent
    return parent_dict

def list_ancestors(parent_dict, node):
    result = []
    while parent_dict[node] != node:
        node = parent_dict[node]
        result.append(node)
    return result

def is_ancestor(parent_dict, ancestor, descendant):
    return ancestor in list_ancestors(parent_dict, descendant)

def classify(parent_dict, groups, taxid, default='Other'):
    for g, name in groups:
        if is_ancestor(parent_dict, g, taxid):
            return name
    return default

#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('taxonomy_file', help='NCBI taxonomy hierarchy file (nodes.dmp) = |-separated table with columns NodeID, ParentID ... ', type=str)
    parser.add_argument('input_file', help='Tab-separated table with columns PDB, TaxID', type=str)
    args = parser.parse_args()
    return vars(args)


def main(taxonomy_file: str, input_file: str) -> None:
    '''Read a tab-separated table domain|taxid from input_file and adds column with superkingdom (Arch/Bact/Euka/Viru), 
    according to NCBI taxonomy hierarchy file taxonomy_file.'''
    with open(taxonomy_file, 'r', encoding=lib.DEFAULT_ENCODING) as f:
        taxonomy = parse_taxonomy_file(f)
    with open(input_file, 'r', encoding=lib.DEFAULT_ENCODING) as f:
        for line in iter(f.readline, ''):
            pdb, taxid, *_ = line.split()
            print(pdb, taxid, classify(taxonomy, GROUPS, taxid), sep='\t')


if __name__ == '__main__':
    args = parse_args()
    main(**args)