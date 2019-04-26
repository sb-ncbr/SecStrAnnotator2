import sys
import json
import argparse

import lib
import no_gap_align

#  CONSTANTS  ##############################################################################

from constants import *

#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('annotations_file', help='JSON file with the list of domains with annotations in SecStrAPI format')
args = parser.parse_args()

annotations_file = args.annotations_file

#  MAIN  ##############################################################################

with open(annotations_file) as r:
    js = json.load(r)

if 'annotations' in js:
    # SecStrAPI format
    domains = [ dom for doms in js['annotations'].values() for name, dom in lib.iterate_names_domains(doms) ]
else:
    # simple format
    domains = list(js.values())

for domain in domains:
    for sse in domain.get(SSES, []):
        label = sse[LABEL]
        magic_number = hash(sse[LABEL]) % 1000
        color =  's' + str(magic_number).zfill(3)
        sse['color'] = color

json.dump(js, sys.stdout, indent=4)
print()
