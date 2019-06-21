import os
from os import path
import sys
import json
import requests
import argparse
from collections import defaultdict
from typing import Tuple

import lib
import no_gap_align

#  CONSTANTS  ##############################################################################

from constants import *

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('all_annotations_file', help='JSON file with the list of domains with annotations in SecStrAPI format')
args = parser.parse_args()

all_annotations_file = args.all_annotations_file

#  MAIN  ##############################################################################

with open(all_annotations_file) as r:
    all_annotations = json.load(r)

pdb2domains = all_annotations[ANNOTATIONS]
for pdb, domains in pdb2domains.items():
    dom_list = domains.values() if isinstance(domains, dict) else domains[:]
    for domain in dom_list:
        for sse in domain[SSES]:
            label = sse.get(LABEL, '')
            metric = sse.get('metric_value')
            print(pdb + domain['chain_id'], label, metric, sep='\t')