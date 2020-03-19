import sys
import json
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--minify', help='Output JSON in minified format (no whitespace) instead of pretty', action='store_true')
args = parser.parse_args()


inp = json.load(sys.stdin)

if args.minify:
    json.dump(inp, sys.stdout, separators=(',', ':'))
else:
    json.dump(inp, sys.stdout, indent=4)