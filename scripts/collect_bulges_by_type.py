import sys
import json
from collections import defaultdict

by_type = defaultdict(list)

for bulge in sys.stdin:
    pdb, typ, chain, start, end = bulge.split()
    by_type[typ].append((pdb, chain, start, end))

sorted_types_bulges = [kv for kv in sorted(by_type.items(), key = lambda kv: len(kv[1]), reverse=True)]

json.dump(dict(sorted_types_bulges), sys.stdout, indent=4)
print()

for typ, bulges in sorted_types_bulges:
    print(typ, len(bulges), file=sys.stderr)
