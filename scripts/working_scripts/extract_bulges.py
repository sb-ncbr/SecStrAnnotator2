import argparse
import os
from os import path
from glob import glob
import json
import sys

import constants as cnst

EXTENSION = '-detected.sses.json'
BULGE_SSE_TYPES = 'bA bP n N m M u U t T s S o O p P q Q r R l L'.split()

parser = argparse.ArgumentParser()
parser.add_argument('in_directory', help='Directory with .sses.json files')
args = parser.parse_args()

in_directory = args.in_directory

files = glob(f'{in_directory}/????{EXTENSION}')
files.sort()

print(len(files), file=sys.stderr)

bulge_types = set()

for i, file in enumerate(files):
    pdb = path.basename(file)[:4]
    with open(file) as f:
        annot = json.load(f)
        sses = annot[pdb][cnst.SSES]
        # print(pdb, len(sses))
    for sse in sses:
        nesteds = sse.get(cnst.NESTED_SSES, ())
        for nested in nesteds:
            # typ = 
            typ = nested[cnst.TYPE]
            if typ in BULGE_SSE_TYPES:
                print(pdb, nested[cnst.COMMENT], nested[cnst.CHAIN_ID], nested[cnst.START], nested[cnst.END])
                bulge_types.add(nested[cnst.COMMENT])
    if i%1000 == 0:
        print(i, file=sys.stderr)

print('Bulge types:', file=sys.stderr)
print(*sorted(bulge_types), file=sys.stderr)