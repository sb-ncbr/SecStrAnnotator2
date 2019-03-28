import json
import os
from os import path
import sys
from pymol import cmd
from pymol import util
cmd.set('cif_use_auth', False)

# How to run: pymol -qcr load_clustered_sses_cif.py -- STRUCTURE_DIRECTORY  CONSENSUS_STRUCTURE CONSENSUS_SSES OUT_SESSION

if len(sys.argv) != 5:
	print('Usage: '+sys.argv[0]+'  STRUCTURE_DIRECTORY  CONSENSUS_STRUCTURE  CONSENSUS_SSES  OUT_SESSION')
	print('STRUCTURE_DIRECTORY and CONSENSUS_STRUCTURE can be empty strings.')
	print(len(sys.argv))
	exit()
_, directory, consensus_structure, consensus_sses, out_session = sys.argv

BASE_COLOR = 'gray80'
DOMAIN = 1 # index of domain name in tuples in samples
MIN_DASH_RADIUS = 0.1  # when radius is dependent on occurrence
MAX_DASH_RADIUS = 1.5
DEFAULT_DASH_RADIUS = 0.25

def load_consensus(consensus_structure, consensus_sses):
	cmd.load(path.join(consensus_structure), 'cons')
	# with open(path.join(directory, 'sample.json')) as f:
	# 	sample = json.load(f)
	# ref_domain = sample[0][DOMAIN]
	# for i, (pdb, domain, chain, rang) in enumerate(sample):
	# 	print(domain)
	# 	cmd.load(path.join(directory, domain+'.cif'))
	# 	color_by_annotation(domain, path.join(directory, domain+'-clust.sses.json'), BASE_COLOR)
	name = 'cons'
	group_name = segment_group_name(name)
	if path.isfile(consensus_sses):
		with open(consensus_sses) as f:
			consensus = json.load(f)	
		cmd.group(group_name)
		for sse in consensus['consensus']['secondary_structure_elements']:
			create_line_segment(sse, group_name)
	cmd.hide()
	cmd.show('cartoon', name)
	util.chainbow(name)  # cmd.color('white', name)
	cmd.show('dashes', group_name)
	cmd.set('dash_gap', 0)
	cmd.dss('cons')
	cmd.zoom('vis')

def load_all(directory, consensus_sses):
	with open(path.join(directory, 'sample.json')) as f:
		sample = json.load(f)
	ref_domain = sample[0][DOMAIN]
	for i, (pdb, domain, chain, rang) in enumerate(sample):
		print(domain)
		cmd.load(path.join(directory, domain+'.cif'), domain)
		color_by_annotation(domain, path.join(directory, domain+'-clust.sses.json'), BASE_COLOR)
	# if path.isfile(consensus_sses):
	# 	with open(consensus_sses) as f:
	# 		consensus = json.load(f)
	# 	cmd.group(segment_group_name('cons'))
	# 	for sse in consensus['consensus']['secondary_structure_elements']:
	# 		create_line_segment(sse, segment_group_name('cons'))
	cmd.hide()
	cmd.show('cartoon')
	cmd.show('dashes')
	cmd.set('dash_gap', 0)
	# cmd.dss()
	cmd.zoom('vis')

def color_by_annotation(domain, annotation_file, base_color, show_line_segments=True):
	with open(annotation_file) as f:
		sses = json.load(f)[domain]['secondary_structure_elements']
	cmd.color(base_color, domain + ' & symbol c')
	cmd.group(sses_group_name(domain))
	if show_line_segments:
		cmd.group(segment_group_name(domain))
	for sse in sses:
		label = sse['label']
		chain = sse['chain_id']
		start = str(sse['start'])
		end = str(sse['end'])
		color = sse['color']
		sel_name = sses_group_name(domain) + '.' + label  # domain + '_' + label
		sel_definition = domain + ' & chain ' + chain + ' & resi ' + start + '-' + end + ' & symbol c'
		if cmd.count_atoms(sel_definition)>0:
			cmd.select(sel_name, sel_definition)
			cmd.color(color, sel_name)
		if show_line_segments:
			create_line_segment(sse, segment_group_name(domain))
		cmd.deselect()
	cmd.dss(domain)

def segment_group_name(domain):
	return domain + '_seg'

def sses_group_name(domain):
	return domain + '_sses'

def create_line_segment_old(sse, domain_name):
	label = sse['label']
	start_vector = sse['start_vector']
	end_vector = sse['end_vector']
	color = sse['color']
	radius = sse['occurrence'] * (MAX_DASH_RADIUS-MIN_DASH_RADIUS) + MIN_DASH_RADIUS if 'occurrence' in sse else DEFAULT_DASH_RADIUS
	cmd.pseudoatom('start', pos=start_vector)
	cmd.pseudoatom('end', pos=end_vector)
	distance_name = domain_name + '_' + label + '_seg'
	cmd.distance(distance_name, 'start', 'end')
	cmd.color(color, distance_name)
	cmd.set('dash_radius', radius, distance_name)
	if sse['type'] in 'GHIh':
		cmd.set('dash_round_ends', 0, distance_name)
	cmd.delete('start')
	cmd.delete('end')

def create_line_segment(sse, group_name):
	label = sse['label']
	start_vector = sse['start_vector']
	end_vector = sse['end_vector']
	color = sse['color']
	radius = sse['occurrence'] * (MAX_DASH_RADIUS-MIN_DASH_RADIUS) + MIN_DASH_RADIUS if 'occurrence' in sse else DEFAULT_DASH_RADIUS
	cmd.pseudoatom('start', pos=start_vector)
	cmd.pseudoatom('end', pos=end_vector)
	distance_name = group_name + '.' + label
	cmd.distance(distance_name, 'start', 'end')
	cmd.color(color, distance_name)
	cmd.set('dash_radius', radius, distance_name)
	if sse['type'] in 'GHIh':
		cmd.set('dash_round_ends', 0, distance_name)
	cmd.delete('start')
	cmd.delete('end')

if consensus_structure != '':
	load_consensus(consensus_structure, consensus_sses)
if directory != '':
	load_all(directory, consensus_sses)
cmd.save(out_session)