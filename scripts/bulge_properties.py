import sys

lines = sys.stdin.readlines()
for line in lines[::2]:
    typ, count = line.split()
    orientation, side0, side1, *_ = typ.split('_')
    if len(side0) > len(side1) or len(side0) == len(side1) and side0 < side1:
        side0, side1 = side1, side0
    print(orientation, side0, side1, len(side0), len(side1), count)
