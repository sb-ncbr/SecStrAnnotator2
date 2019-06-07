import os
import sys
import json

if len(sys.argv) < 3:
	print('Usage: '+sys.argv[0]+' TAXFILE INPUTFILE')
	print(' TAXFILE = NCBI taxonomy hierarchy file (nodes.dmp) = |-separated table with columns NodeID, ParentID ... ')
	print(' INPUTFILE = Space-separated table with columns PDB, TaxID')
	exit()
taxonomy_file = sys.argv[1]
input_file = sys.argv[2]

EUKARYOTA = '2759'
BACTERIA = '2'
ARCHEA = '2157'
HUMAN = '9606'
VIRUSES = '10239'
GROUPS = [(EUKARYOTA,'Euka'), (BACTERIA,'Bact'), (ARCHEA,'Arch'), (VIRUSES, 'Viru')]

# PHYLA = '976 1090 1117 1224 1239 1297 2833 2836 2870 3035 3041 4761 4890 5204 5747 5794 6029 6040 6073 6157 6217 6231 6340 6447 6656 7568 7586 7711 10190 10197 10205 10214 10219 10226 10229 10232 27563 28889 28890 31291 32066 33209 33310 33313 33467 35493 40117 42241 43120 51516 51967 57723 62680 65842 66780 67810 67812 67814 67817 67818 67819 67820 68297 69815 74015 74152 74201 91989 95818 104731 134625 142182 142187 157124 192989 200783 200795 200918 200930 200938 200940 201174 203682 203691 204428 221216 256845 265317 310840 363464 419944 422282 447827 451459 456828 508458 544448 569577 640293 651137 743724 743725 877183 928852 1031332 1052197 1052815 1104542 1134404 1154675 1154676 1293497 1312402 1379697 1383058 1448051 1448933 1448937 1462422 1462430 1492816 1618330 1618338 1618339 1618340 1619053 1655434 1696033 1703755 1704031 1706441 1729712 1752708 1752716 1752717 1752718 1752719 1752720 1752721 1752722 1752723 1752724 1752725 1752726 1752727 1752728 1752729 1752730 1752731 1752732 1752733 1752734 1752735 1752736 1752737 1752738 1752739 1752740 1752741 1798710 1801616 1801631 1802339 1817796 1817797 1817798 1817799 1817800 1817801 1817802 1817803 1817804 1817805 1817806 1817807 1817808 1817809 1817810 1817811 1817812 1817898 1817899 1817900 1817901 1817902 1817903 1817904 1817905 1817906 1817907 1817908 1817909 1817910 1817911 1817912 1817913 1817914 1817915 1817916 1817917 1817918 1817919 1817920 1817921 1817922 1817923 1817924 1817925 1817926 1819803 1853220 1855361 1913637 1913638 1915410 1930617 1936271 1936272 1936987 2013583'
# PHYLA = PHYLA.split()
# phylum_names = dict(zip(PHYLA, PHYLA))
# phylum_names.update({ '1117': 'Cyanobacteria', '1224': 'Proteobacteria', '1239': 'Firmicutes', '1297': 'Deinococcus-Thermus', '201174': 'Actinobacteria', 
# 	'28889': 'Crenarchaeota', '28890': 'Euryarchaeota', '35493': 'Streptophyta', '4890': 'Ascomycota', '7711': 'Chordata', 'Other': 'Other' })
# GROUPS = list(phylum_names.items())


def read_json(filename):
	with open(filename, 'r') as f:
		result = json.load(f)
	return result

def parse_taxonomy_file(f):
	parent_dict = {}
	for line in iter(f.readline, ''):
		node, parent, *_ = (x.strip() for x in line.split('|'))
		parent_dict[node] = parent
	return parent_dict

def list_ancestors(parent_dict, node):
	result = []
	while parent_dict[node] != node:
		node = parent_dict[node]
		result.append(node)
	return result

def is_ancestor(parent_dict, ancestor, descendant):
	return ancestor in list_ancestors(parent_dict, descendant)

def classify(parent_dict, groups, taxid, default='Other'):
	for g, name in groups:
		if is_ancestor(parent_dict, g, taxid):
			return name
	return default

with open(taxonomy_file) as f:
	taxonomy = parse_taxonomy_file(f)

with open(input_file) as f:
	for line in iter(f.readline, ''):
		pdb, taxid, *_ = line.split()
		#print(pdb, taxid, ' '.join(list_ancestors(taxonomy, taxid)))
		print(pdb, taxid, classify(taxonomy, GROUPS, taxid), sep='\t')

exit()


if nonred:
	pdblist = sorted(sorted(uni2pdb[uni])[0] for uni in uni2pdb)
else:
	pdblist = sorted(pdb for uni in uni2pdb for pdb in uni2pdb[uni])

if os.path.isfile(outdir):
	raise Exception('Cannot create output directory '+outdir+' because a file with this name exists.')
if not os.path.isdir(outdir):
	os.makedirs(outdir)
os.chdir(outdir)

print('Fetching ' + ('nonredundant' if nonred else 'all') + ' PDBs:')
print(str(len(pdblist)) + ' PDBs')

command = 'fetch ' + ' '.join(pdblist) + ', type=pdb, async=0'

os.system('pymol -Qcd \'' + command + '\'')

failed = []
for pdb in pdblist:
	if not os.path.isfile(pdb+'.pdb'):
		failed.append(pdb)
print('Failed to fetch ' + str(len(failed)) + ' PDBs: ' + ', '.join(failed))


