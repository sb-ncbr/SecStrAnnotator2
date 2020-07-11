import os
import sys
import json
import requests
import argparse
from collections import defaultdict

import lib

#  CONSTANTS  ##############################################################################

from constants import *


#  PARSE ARGUMENTS  ##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('publication_info', help='JSON file from get_publication_info.py', type=str)
parser.add_argument('--condensed', help='Merge result by article name', action='store_true')
parser.add_argument('--primary', help='Take only the primary citation (first listed) for each PDB', action='store_true')
parser.add_argument('--domain_list', help='JSON file with the list of domains in SecStrAPI format (from merge_domain_lists.py) to include UniProt names in output', type=str, default=None)
args = parser.parse_args()

publinfo_file = args.publication_info
condensed = args.condensed
primary = args.primary
domain_list_file = args.domain_list

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

if domain_list_file is not None:
	with open(domain_list_file, 'r') as r:
		domain_list = json.load(r)
	pdb2uni = {}
	for pdb, doms in domain_list['annotations'].items():
		uni = next(iter(doms.values()))['uniprot_name']
		pdb2uni[pdb] = uni
else:
	pdb2uni = None

by_title = defaultdict(lambda: [])

columns = ['PDB', 'Year', 'Author', 'Title', 'DOI']
if pdb2uni is not None:
	columns.append('UniProtName')
print(*columns, sep='\t')

for pdb, pubs in publinfo.items():
	if primary:
		pubs = pubs[:1]
	for pub in pubs:
		title = pub['title'].replace('\n', ' ').replace('\t', ' ')
		year = pub['journal_info']['year']
		author = pub['author_list'][0]['last_name']
		if author is None:
			author =  pub['author_list'][0]['full_name']
		doi = pub['doi']
		if not condensed:
			outs = [pdb, year, author, title, doi]
			if pdb2uni is not None:
				outs.append(pdb2uni[pdb])
			print(*outs, sep='\t')
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
			doi = dois[0]
			# raise Exception(f'Ambiguous DOI: {dois}')
		if len(years) == 0:
			year = None
		elif len(years) == 1:
			year = years[0]
		else:
			raise Exception(f'Ambiguous year: {years}')
		title = titles[0]
		author = min(authors, key=len)
		outs = [','.join(pdbs), year, author, title, doi]
		if pdb2uni is not None:
			unis = set( pdb2uni[pdb] for pdb in pdbs )
			unis.discard(None)
			unis = sorted(unis)
			if len(unis) !=1:
				sys.stderr.write(f'More UniProts in one paper: {pdbs} -> {unis}\n')
			outs.append(','.join(unis))				
		print(*outs, sep='\t')
