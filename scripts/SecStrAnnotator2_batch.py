import json
import os
from os import path
import argparse
import sys
from collections import OrderedDict
try: 
	import queue
except ImportError:
	import Queue as queue
import threading
import subprocess

#  PARSE ARGUMENTS  ###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('directory', help='directory with the input and output files (argument DIRECTORY for SecStrAnnotator)', type=str)
parser.add_argument('template', help='template domain specification (argument TEMPLATE for SecStrAnnotator)', type=str)
parser.add_argument('queries', help='JSON file with the list of domains to be annotated (in format {PDB:[[domain_name,chain,ranges]]}, will be processed to QUERY arguments for SecStrAnnotator)', type=str)
parser.add_argument('--options', help='Any options that are to be passed to SecStrAnnotator (must be enclosed in quotes and contain spaces, not to be confused with Python arguments, e.g. --options \'--ssa dssp --soft\' or \' --soft\')', type=str, default='')
parser.add_argument('--threads', help='Number of parallel threads (default: 1)', type=int, default=1)
parser.add_argument('--dll', help='Path to the SecStrAnnotator DLL (default: SecStrAnnotator2.dll)', type=str, default='SecStrAnnotator2.dll')
args = parser.parse_args()

options = [opt for opt in args.options.split(' ') if opt != '']
onlyssa = '--onlyssa' in options
directory = args.directory
template = args.template
queries_file = args.queries
n_threads = args.threads
RUN_DIR = path.dirname(__file__)  # os.getcwd()
SECSTRANNOTATOR = path.join(RUN_DIR, args.dll)

all_annotations_file = path.join(directory, 'all_annotations.sses.json')
output = path.join(directory, 'stdout.txt')
output_err = path.join(directory, 'stderr.txt')
out_files_extensions = ['-aligned.cif', '-detected.sses.json', '-annotated.sses.json', '-annotated.pse']

#  FUNCTIONS  ###############################################################################

def clear_file(filename):
	with open(filename, 'w') as w:
		w.write('')

def try_rename_file(source, dest):
	try:
		os.rename(source, dest)
	except:
		pass

def copy_file(source, dest, append=False):
	with open(source, 'r') as r:
		with open(dest, 'a' if append else 'w') as w:
			for line in iter(r.readline, ''):
				w.write(line)

def remove_file(filename):
	os.remove(filename)

def try_read_json(filename):
	with open(filename) as f:
		try:
			result = json.load(f, object_pairs_hook=OrderedDict)
		except ValueError as e:
			sys.stderr.write('Error: "' + filename + '" is not a valid JSON file (' + str(e) +  ') \n')
			exit(1)
	return result

def check_json_type(json_object, typeex):
	if isinstance(typeex, tuple):
		return any( check_json_type(json_object, typ) for typ in typeex )
	elif isinstance(typeex, str):
		return isinstance(json_object, str)
	elif isinstance(typeex, bool):
		return isinstance(json_object, bool)
	elif isinstance(json_object, (int, float)):
		return isinstance(json_object, (int, float))
	elif typeex is None:
		return json_object is None
	elif isinstance(typeex, list):
		if not isinstance(json_object, list):
			return False 
		if len(typeex) > 0:
			return all( check_json_type(elem, typeex[0]) for elem in json_object )
		else: 
			return True
	elif isinstance(typeex, dict):
		if not isinstance(json_object, dict):
			return False 
		if len(typeex) > 0:
			value_typeex = next(iter(typeex.values()))
			return all( isinstance(key, str) and check_json_type(value, value_typeex) for key, value in json_object.items() )
		else: 
			return True
	raise Exception('Invalid typeex')

def read_domains_json(filename):
	result = try_read_json(filename)
	if check_json_type(result, {'PDB':[['']]}):
		return result
	else:
		sys.stderr.write('Error: Expected a JSON format {PDB:[[domain_name,chain,ranges]]} in "' + filename + '" \n')
		exit(1)

class ProgressBar:
	def __init__(self, n_steps, width=100, title='', writer=sys.stdout):
		self.n_steps = n_steps # expected number of steps
		self.width = width
		self.title = (' '+title+' ')[0:min(len(title)+2, width)]
		self.writer = writer
		self.done = 0 # number of completed steps
		self.shown = 0 # number of shown symbols
	def start(self):
		self.writer.write('|' + self.title + '_'*(self.width-len(self.title)) + '|\n')
		self.writer.write('|')
		self.writer.flush()
		return self
	def step(self, n_steps=1):
		if self.n_steps == 0:
			return
		self.done = min(self.done + n_steps, self.n_steps)
		new_shown = int(self.width * self.done / self.n_steps)
		self.writer.write('*' * (new_shown-self.shown))
		self.writer.flush()
		self.shown = new_shown
	def finalize(self):
		self.step(self.n_steps - self.done)
		self.writer.write('|\n')
		self.writer.flush()

def run_in_threads (do_job, jobs, n_threads, progress_bar=None, initialize_thread_sync=None, finalize_thread_sync=None):
	q = queue.Queue()
	for job in jobs:
		q.put(job)
	def worker():
		all_done = False
		while not all_done:
			try:
				job = q.get(block=False)
				do_job(job)
				if progress_bar is not None:
					progress_bar.step()
				q.task_done()
			except queue.Empty:
				all_done = True
	# threads = []
	# for i in range(n_threads):
	# 	thread = threading.Thread(target=worker)
	# 	threads.append(thread)
	# 	if initialize_thread_sync is not None:
	# 		initialize_thread_sync(thread)
	threads = [ threading.Thread(target=worker) for i in range(n_threads) ]
	if progress_bar is not None:
		progress_bar.start()
	if initialize_thread_sync is not None:
		for thread in threads:
			initialize_thread_sync(thread)
	for thread in threads:
		thread.start()
	q.join()
	if finalize_thread_sync is not None:
		for thread in threads:
			finalize_thread_sync(thread)
	if progress_bar is not None:
		progress_bar.finalize()
	
#  MAIN  ###############################################################################

# Determine whether can run dotnet SecStrAnnotator.dll and whether the template files exist
if not path.isfile(SECSTRANNOTATOR):
	sys.stderr.write('Error: "' + SECSTRANNOTATOR + '" not found\n')
	exit(1)

template_pdb = template.split(',')[0]
template_struct_file = path.join(directory, template_pdb+'.cif')
template_annot_file = path.join(directory, template_pdb+'-template.sses.json')
if not path.isfile(template_struct_file):
	sys.stderr.write('Error: "' + template_struct_file + '" not found\n')
	exit(1)
if not path.isfile(template_annot_file):
	sys.stderr.write('Error: "' + template_annot_file + '" not found\n')
	exit(1)

secstrannotator_commands = ['dotnet', SECSTRANNOTATOR]

# Prepare domain list
domains = read_domains_json(queries_file)
pdbs = sorted(domains)
if not path.isdir(directory):
	sys.stderr.write('Error: "' + directory + '" is not a directory\n')
	exit(1)
print(directory)
found_pdbs = [pdb for pdb in pdbs if path.isfile(path.join(directory, pdb+'.cif'))]
not_found_pdbs = [pdb for pdb in pdbs if not path.isfile(path.join(directory, pdb+'.cif'))]
n_domains = sum(len(domains[pdb]) for pdb in pdbs)
n_found_domains = sum(len(domains[pdb]) for pdb in found_pdbs)
not_found_domains = [domain for pdb in not_found_pdbs for domain in domains[pdb] ]

print('Listed ' + str(len(pdbs)) + ' PDBs (' + str(n_domains) + ' domains), found ' + str(len(found_pdbs)) + ' PDBs (' + str(n_found_domains) + ' domains)')

all_annotations = OrderedDict()
failed = []

def process_pdb(pdb):
	thread_out = path.join(directory, 'stdout_thread_' + threading.current_thread().name + '.txt')
	thread_err = path.join(directory, 'stderr_thread_' + threading.current_thread().name + '.txt')
	clear_file(path.join(directory, pdb + '.label2auth.tsv'))
	for domain_name, chain, ranges in domains[pdb]:
		query = pdb + ',' + chain + ',' + ranges
		with open(thread_out, 'a') as w_out:
			with open(thread_err, 'a') as w_err:
				w_out.write('\n' + '-'*70 + '\n' + domain_name + '\n')
				w_err.write('-'*70 + '\n' + domain_name + '\n')
				w_out.flush()
				w_err.flush()
				regular_arguments = [directory, template, query] if not onlyssa else [directory, query]
				exit_code = subprocess.call(secstrannotator_commands + regular_arguments + options, stdout=w_out, stderr=w_err)
		if exit_code == 0:
			for ext in out_files_extensions:
				try_rename_file(path.join(directory, pdb + ext), path.join(directory, domain_name + ext))
			try:
				copy_file(path.join(directory, pdb + '-label2auth.tsv'), path.join(directory, pdb + '.label2auth.tsv'), append=True)
				remove_file(path.join(directory, pdb + '-label2auth.tsv'))
			except:
				pass
			if not onlyssa:
				annotation = try_read_json(path.join(directory, domain_name + '-annotated.sses.json')).get(pdb, {})
				annotation['domain'] = { 'name': domain_name, 'pdb': pdb, 'chain': chain, 'ranges': ranges }
				all_annotations[domain_name] = annotation
		else:
			failed.append(domain_name)

def clear_outputs(thread):
	thread_out = path.join(directory, 'stdout_thread_' + thread.name + '.txt')
	thread_err = path.join(directory, 'stderr_thread_' + thread.name + '.txt')
	clear_file(thread_out)
	clear_file(thread_err)

def merge_outputs(thread):
	thread_out = path.join(directory, 'stdout_thread_' + thread.name + '.txt')
	thread_err = path.join(directory, 'stderr_thread_' + thread.name + '.txt')
	copy_file(thread_out, output, append=True)
	copy_file(thread_err, output_err, append=True)
	remove_file(thread_out)
	remove_file(thread_err)

clear_file(output)
clear_file(output_err)

progress_bar = ProgressBar(len(found_pdbs), title='Running SecStrAnnotator on '+str(n_found_domains)+' domains', width=70)

# Run SecStrAnnotator in parallel threads
run_in_threads(process_pdb, found_pdbs, n_threads, progress_bar=progress_bar, initialize_thread_sync=clear_outputs, finalize_thread_sync=merge_outputs)

# Output collected data
with open(all_annotations_file, 'w') as w:
	json.dump(all_annotations, w, indent=4)

print('Failed to find ' + str(len(not_found_domains)) + ' domains:')
print(', '.join(name for name, chain, ranges in not_found_domains))
print('Failed to annotate ' + str(len(failed)) + ' domains:')
print(', '.join(failed))
