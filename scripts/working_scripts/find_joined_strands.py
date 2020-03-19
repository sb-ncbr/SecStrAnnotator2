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
parser.add_argument('annotation_file', help='JSON file with the list of domains with annotations in SecStrAPI format', type=str)
args = parser.parse_args()

annotation_file = args.annotation_file

#  MAIN  ##############################################################################

def is_sheet(sse):
    return sse['type'] in 'EBe'

with open(annotation_file) as r:
    annotations = json.load(r)

pdb2domains = annotations[ANNOTATIONS]
for pdb, domains in pdb2domains.items():
    for domain, annot in domains.items():
        sses = annot[SSES]
        for sse in sses:
            if is_sheet(sse) and any( is_sheet(nested) for nested in sse.get('nested_sses', []) ):
                print(domain, sse[LABEL])
            elif sse[LABEL] == '4-2':
                if any( sse[LABEL] == '4-1' for sse in sses ) and any( sse[LABEL] == '4-3' for sse in sses ):
                    print(f'4-2 is naturally joined in {pdb}')
            elif sse[LABEL] == '3-2':
                if any( sse[LABEL] == '3-1' for sse in sses ) and any( sse[LABEL] == '3-3' for sse in sses ):
                    print(f'3-2 is naturally joined in {pdb}')