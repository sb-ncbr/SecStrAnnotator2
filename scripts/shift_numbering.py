import json
from collections import OrderedDict



INPUT_FILE = '1og2-template-auth.sses.json'
OUTPUT_FILE = '1og2-template.sses.json'
SHIFT = -19

with open(INPUT_FILE) as f:
    annot = json.load(f, object_pairs_hook=OrderedDict)

for sse in annot['1og2']['secondary_structure_elements']:
    sse['start'] += SHIFT
    sse['end'] += SHIFT

with open(OUTPUT_FILE, 'w') as w:
    json.dump(annot, w)