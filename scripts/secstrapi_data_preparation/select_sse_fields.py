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

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('input_annotations_file', help='JSON file with the list of domains with annotations in SecStrAPI format', type=str)
parser.add_argument('--fields', help='Comma-separated list of fields which should be kept in "secondary_structure_elements"', 
    type=str,
    default='label,chain_id,start,end,auth_chain_id,auth_start,auth_start_ins_code,auth_end,auth_end_ins_code,type,color,metric_value,confidence,sequence'
)
args = parser.parse_args()

input_annotations_file = args.input_annotations_file
fields = args.fields.split(',')

#  MAIN  ##############################################################################

with open(input_annotations_file) as r:
    all_annotations = json.load(r)

pdb2domains = all_annotations[ANNOTATIONS]
for pdb, domains in pdb2domains.items():
    for domain, annot in domains.items():
        for sse in annot[SSES]:
            for field in list(sse.keys()):
                if field not in fields:
                    sse.pop(field)

json.dump(all_annotations, sys.stdout, indent=4)
print()
