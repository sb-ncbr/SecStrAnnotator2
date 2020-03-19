import json
import os
import os.path
import sys
import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('all_annotations_file', help='JSON file with annotations in SecStrAPI format', type=str)
# parser.add_argument('--fields', help='SSE fields to extract as columns', type=str, default='')
parser.add_argument('--one_domain_per_pdb', help='Take only first domain for each PDB', action='store_true')
args = parser.parse_args()

all_annotations_file = args.all_annotations_file
# fields = [ field for field in args.fields.split(',') if field != '' ]
one_domain_per_pdb = args.one_domain_per_pdb


SSES='secondary_structure_elements'
LABEL='label'
TYPE='type'
START='start'
END='end'
CHAIN='chain_id'
COMMENT='comment'
SEQ='sequence'
NESTED='nested_sses'

bulge_types='nNmMtTsSoOpPqQrRlL'

def equal(sse1, sse2):
	return difference(sse1, sse2) == 0
	if sse1 == None or sse2 == None:
		return False
	else:
		return sse1[START] == sse2[START] and sse1[END] == sse2[END]

def difference(sse1, sse2):
	if sse1 == None or sse2 == None:
		return math.nan
	else:
		return abs(sse1[START] - sse2[START]) + abs(sse1[END] - sse2[END])

def length(sse):
	return sse[END] - sse[START] + 1

def longest_nested(sse, nested_type):
	if sse is not None and NESTED in sse:
		nested = sse[NESTED]
		lengths = [length(nes) for nes in nested if nes[TYPE]==nested_type]
		# sys.stderr.write('Lengths: ' + str(lengths) + '\n')
		longest = max(lengths + [0])
		return longest
	else:
		return 0

def longest_nested_starting(sse, nested_type, from_start_tolerance=1):
	if sse is not None and NESTED in sse:
		nested = sse[NESTED]
		lengths = [length(nes) for nes in nested if nes[TYPE]==nested_type and nes[START]<=sse[START]+from_start_tolerance]
		# sys.stderr.write('Lengths: ' + str(lengths) + '\n')
		longest = max(lengths + [0])
		return longest
	else:
		return 0

def longest_nested_ending(sse, nested_type, from_end_tolerance=1):
	if sse is not None and NESTED in sse:
		nested = sse[NESTED]
		lengths = [length(nes) for nes in nested if nes[TYPE]==nested_type and nes[END]>=sse[END]-from_end_tolerance]
		# sys.stderr.write('Lengths: ' + str(lengths) + '\n')
		longest = max(lengths + [0])
		return longest
	else:
		return 0

def elemwise_sum(lists):
	return [ sum(elems) for elems in zip(*lists) ]

def count_GHI_bonds(sse):
	if sse is None:
		return [0, 0, 0]
	elif NESTED in sse:
		return elemwise_sum( count_GHI_bonds(nested) for nested in sse[NESTED] )
	else:
		if sse[TYPE] == 'G':
			return [length(sse) - 1, 0, 0]
		elif sse[TYPE] == 'H':
			return [0, length(sse) - 2, 0]
		elif sse[TYPE] == 'I':
			return [0, 0, length(sse) - 3]
		else:
			return [0, 0, 0]

def sse_to_row(sse, *extras):
	row = []
	row.extend( extra if extra is not None else 'NA' for extra in extras )
	row.extend( sse[field] if sse is not None else 'NA' for field in fields )
	row.extend([longest_nested(sse, 'G'), longest_nested_starting(sse, 'G'), longest_nested_ending(sse, 'G'), 
		longest_nested(sse, 'H'), longest_nested(sse, 'I'), longest_nested_starting(sse, 'I'), longest_nested_ending(sse, 'I'), *count_GHI_bonds(sse)])
	return row

def find_bulges_recursive(sse):
	result = []
	if sse[TYPE] in bulge_types:
		result.append(sse)
	for nested in sse.get(NESTED, []):
		result.extend(find_bulges_recursive(nested))
	return result
	#TODO write as generator


with open(all_annotations_file, 'r') as f:
	annotations=json.load(f)['annotations']

labels = sorted(set( sse[LABEL] for pdb_annot in annotations.values() for dom_annot in pdb_annot.values() for sse in dom_annot[SSES] ))
# print(labels)
# exit()
result = []

for pdb, pdb_annot in annotations.items():
	domain_annots = list(pdb_annot.items())
	if one_domain_per_pdb:
		domain_annots = domain_annots[:1]
	for domain, domain_annot in domain_annots:
		uni = domain_annot['uniprot_id']
		sse_bulges = [ (sse, bulge) for sse in domain_annot[SSES] for bulge in find_bulges_recursive(sse) ]
		for sse, bulge in sse_bulges:
			row = [uni, pdb, domain, sse[LABEL], bulge[LABEL], bulge[CHAIN], bulge[START], bulge[END], length(bulge), bulge[TYPE]]
			# for field in fields:
			# 	row.append(bulge[field] if bulge!=None and field in bulge else 'NA')
			result.append(row)
				
print('UniProt', 'PDB', 'Domain', 'label', 'bulge_label', CHAIN, 'start', 'end', 'length', 'type', sep='\t')
for row in sorted(result):
	print(*row, sep='\t')















