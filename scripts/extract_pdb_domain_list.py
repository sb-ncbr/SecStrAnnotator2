import sys
import json
import argparse

import lib

#  CONSTANTS  ##############################################################################

from constants import *

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='Domain list in SecStrAPI format', type=str)
args = parser.parse_args()

input_file = args.input_file

#  MAIN  ##############################################################################

with open(input_file) as r:
    domain_list = json.load(r)
    pdb2domains = domain_list[ANNOTATIONS]

simple_list = {}
for pdb, domains in pdb2domains.items():
    simple_list[pdb] = [ name for name, dom in lib.iterate_names_domains(domains) ] 

json.dump(simple_list, sys.stdout, indent=4)
print()

n_pdbs = len(simple_list)
n_domains = sum( len(doms) for doms in simple_list.values() )
sys.stderr.write(f'Extracted a list of {n_domains} domains in {n_pdbs} PDB entries\n')