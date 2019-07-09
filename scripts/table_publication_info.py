import os
import sys
import json
import requests
import argparse
from collections import defaultdict

#  CONSTANTS  ##############################################################################

from constants import *


#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('publication_info', help='JSON file from get_publication_info.py', type=str)
parser.add_argument('--condensed', help='Merge result by article name', action='store_true')
args = parser.parse_args()

publinfo_file = args.publication_info
condensed = args.condensed

#  MAIN  ##############################################################################

def unique_nonnull_zip(*tuples):
	k = len(tuples[0])
	results = tuple( set() for i in range(k) )
	for t in tuples:
		for i, x in enumerate(t):
			if x is not None:
				results[i].add(x)
	return tuple( sorted(res) for res in results )

def only_alnum(string):
	return ''.join( c for c in string if c.isalnum() )

with open(publinfo_file, 'r') as r:
	publinfo = json.load(r)

by_title = defaultdict(lambda: [])

print('PDB', 'Year', 'Author', 'Title', 'DOI', sep='\t')

for pdb, pubs in publinfo.items():
	for pub in pubs:
		title = pub['title'].replace('\n', ' ').replace('\t', ' ')
		year = pub['journal_info']['year']
		author = pub['author_list'][0]['last_name']
		if author is None:
			author =  pub['author_list'][0]['full_name']
		doi = pub['doi']
		if not condensed:
			print(pdb, year, author, title, doi, sep='\t')
		by_title[only_alnum(title).upper()].append((pdb, year, author, title, doi))

if condensed:
	for title, pubs in by_title.items():
		pdbs, years, authors, titles, dois = unique_nonnull_zip(*pubs)
		if len(dois) == 0:
			doi = None
		elif len(dois) == 1:
			doi = dois[0]
		else:
			sys.stderr.write(f'Ambiguous DOI: {dois}\n')
			doi = doi[0]
			# raise Exception(f'Ambiguous DOI: {dois}')
		if len(years) == 0:
			year = None
		elif len(years) == 1:
			year = years[0]
		else:
			raise Exception(f'Ambiguous year: {years}')
		title = titles[0]
		author = min(authors, key=len)
		print( ','.join(pdbs), year, author, title, doi, sep='\t')
