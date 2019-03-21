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

OUTPUT_EXT = '.json'

#  FUNCTIONS  ##############################################################################

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('all_annotations_file', help='JSON file with the list of domains with annotations in SecStrAPI format', type=str)
parser.add_argument('output_directory', help='Directory to output per-PDB-entry annotations', type=str)
args = parser.parse_args()

all_annotations_file = args.all_annotations_file
output_directory = args.output_directory

#  MAIN  ##############################################################################

with open(all_annotations_file) as r:
    all_annotations = json.load(r)

shutil.rmtree(output_directory, ignore_errors=True)
os.makedirs(output_directory)

api_version = all_annotations[API_VERSION]
pdb2domains = all_annotations[ANNOTATIONS]
for pdb, domains in pdb2domains.items():
    result = { API_VERSION: api_version, ANNOTATIONS: { pdb: domains } }
    with open(path.join(output_directory, pdb + OUTPUT_EXT), 'w') as w:
        json.dump(result, w, indent=4)

n_pdbs = len(all_annotations[ANNOTATIONS])
n_domains = sum( len(doms) for doms in all_annotations[ANNOTATIONS].values() )
sys.stderr.write(f'Collected annotations for {n_domains} domains in {n_pdbs} PDB entries\n')