import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument('in_file')
parser.add_argument('out_file')
parser.add_argument('--times', type=float)
args = parser.parse_args()

with open(args.in_file) as r:
    lines = r.readlines()

rand = random.Random()
line_indices = range(len(lines))
for i in range(int(args.times * len(lines))):
    x, y = rand.choices(line_indices, k=2)
    lines[x], lines[y] = lines[y], lines[x]

with open(args.out_file, 'w') as w:
    w.writelines(lines)