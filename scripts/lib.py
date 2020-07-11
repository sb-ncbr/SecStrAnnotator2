import sys
import json
import requests
import shutil
import ftplib
import tarfile
from collections import namedtuple
from typing import List, Dict, Tuple, Union, Iterator, Iterable, Callable, TypeVar, Optional

from constants import *

DEFAULT_ENCODING = 'utf-8'

K = TypeVar('K')
V = TypeVar('V')

DEFAULT_STRUCTURE_QUALITY = 0.0

def insert_after(dictionary: Dict[K, V], after_what: K, new_key_value_pairs: Iterable[Tuple[K, V]]) -> None:
    key_value_pairs = list(dictionary.items())
    dictionary.clear()
    for key, value in key_value_pairs:
        dictionary[key] = value
        if key == after_what:
            for k, v in new_key_value_pairs:
                dictionary[k] = v

def iterate_names_domains(domains: Union[Dict[str, V], List[V]]) -> Iterator[Tuple[str, V]]:
    if isinstance(domains, dict):
        yield from domains.items()
    elif isinstance(domains, list):
        for dom in domains:
            yield f'{dom[PDB]},{dom[CHAIN]},{dom[RANGES]}', dom
    else:
        raise Exception('Domains must by stored in a dict or list')
    # name_domain_gen = domains.items() if isinstance(domains, dict) else ( (f'{dom[PDB]},{dom[CHAIN]},{dom[RANGES]}', dom) for dom in domains )

# def new_mapping(domain_name, source, family, pdb, chain, ranges):
#     mapping = { DOMAIN_NAME: domain_name, SOURCE: source, FAMILY_ID: family, PDB: pdb, CHAIN: chain, RANGES: ranges }
#     return mapping

def create_domain_id(pdb: str, chain: str) -> str:
    return pdb + chain

def get_structure_quality(pdb: str) -> float:
    QUALITY_API_URL = 'https://www.ebi.ac.uk/pdbe/api/validation/summary_quality_scores/entry/{pdb}'
    QUALITY = 'overall_quality'
    url = QUALITY_API_URL.format(pdb=pdb)
    r = json.loads(requests.get(url).text)
    quality = r.get(pdb, {}).get(QUALITY, DEFAULT_STRUCTURE_QUALITY)
    return quality

def single(iterable: Iterable[V]) -> V:
    '''Return the single element of the iterable, or raise ValueError if len(iterable) is not 1'''
    iterator = iter(iterable)
    try:
        value = next(iterator)
        try:
            _ = next(iterator)
            raise ValueError('Iterable contains more than one element)')
        except StopIteration:
            return value
    except StopIteration:
        raise ValueError('Iterable contains no elements')

def get_from_ftp(address: str, destination_file: str, username: str = 'anonymous') -> None:
    if address.startswith('ftp://'):
        address = address[6:]
    server, filename = address.split('/', maxsplit=1)
    with ftplib.FTP(server, username) as connection:
        with open(destination_file, 'wb') as w:
            connection.retrbinary(f'RETR {filename}', w.write)

def extract_from_tar(tar_archive: str, extracted_file: str, destination_file: str) -> None:
    with tarfile.open(tar_archive, 'r') as tar:
        with tar.extractfile(extracted_file) as r:
            with open(destination_file, 'wb') as w:
                shutil.copyfileobj(r, w)

class LazyDict:
    def __init__(self, initializer: Callable[[K], V]):
        self.initializer = initializer
        self.dictionary = {}
    def get(self, key: K) -> V:
        if key not in self.dictionary:
            self.dictionary[key] = self.initializer(key)
        return self.dictionary[key]
    def __getitem__(self, key: K) -> V:
        return self.get(key)


class Label2AuthConverter:
    def __init__(self, conversion_table_file: str, unknown_ins_code_as_empty_string: bool=False):
        self.file = conversion_table_file
        self.table = {}
        self.chain_table = {}
        self.unknown_ins_code_as_empty_string = unknown_ins_code_as_empty_string
        with open(conversion_table_file, 'r', encoding=DEFAULT_ENCODING) as r:
            for line in r:
                if not line.lstrip().startswith('#'):
                    chain, resi, auth_chain, auth_resi, auth_ins_code, label_comp_id = line.strip().split('\t')
                    if self.unknown_ins_code_as_empty_string and auth_ins_code == '?':
                        auth_ins_code = ''
                    self.table[(chain, int(resi))] = (auth_chain, int(auth_resi), auth_ins_code)
                    if chain not in self.chain_table:
                        self.chain_table[chain] = auth_chain
                    elif self.chain_table[chain] != auth_chain:
                        raise Exception(f'label_asym_id {chain} maps to multiple auth_asym_id ({self.chain_table[chain]}, {auth_chain})')
        # if '1cpt' in conversion_table_file:
        #     sys.stderr.write(f'CHAIN_TABLE: {self.chain_table}')
        #     sys.stderr.write(f'TABLE: {self.table}')

    def auth_chain(self, chain: str) -> str:
        try:
            auth_chain = self.chain_table[chain]
            return auth_chain
        except KeyError:
            raise Exception(f'Did not find chain {chain} in {self.file}')

    def auth_chain_resi_ins(self, chain: str, resi: int) -> Tuple[str, int, str]:
        # sys.stderr.write(f'{self.table}\n')
        try:
            return self.table[(chain, resi)]
        except KeyError:
            raise KeyError(f'Did not find residue {chain} {resi} in {self.file}')

    def auth_chain_ranges(self, chain: str, ranges: str) -> Tuple[str, str]:
        auth_chain = self.auth_chain(chain)
        auth_ranges = ''
        for rang in ranges.split(','):
            fro, to = rang.split(':')
            if fro == '':
                auth_fro = '' 
            else:
                c, resi, ins = self.auth_chain_resi_ins(chain, int(fro))
                auth_fro = str(resi) + ins
            if to == '':
                auth_to = '' 
            else:
                c, resi, ins = self.auth_chain_resi_ins(chain, int(to))
                auth_to = str(resi) + ins
            auth_ranges += auth_fro + ':' + auth_to
        return auth_chain, auth_ranges


class UniProtManager:
    UNIPROT = 'UniProt'
    UNIPROT_NAME = 'name'
    MAPPINGS = 'mappings'
    CHAIN = 'struct_asym_id'
    def __init__(self, source_url='http://www.ebi.ac.uk/pdbe/api/mappings/uniprot/'):
        self.pdbchain2uniidname = {}
        self.source_url = source_url
    def get_uniprot_id_and_name(self, pdb: str, chain: str) -> Tuple[str, str]:  
        if (pdb, chain) not in self.pdbchain2uniidname:
            self._add_pdb(pdb)
            if (pdb, chain) not in self.pdbchain2uniidname:
                self.pdbchain2uniidname[(pdb, chain)] = (None, None)
        return self.pdbchain2uniidname[(pdb, chain)]
    def _add_pdb(self, pdb: str) -> None:
        url = f'{self.source_url}/{pdb}'
        r = json.loads(requests.get(url).text)
        unip = r[pdb][self.UNIPROT]
        for uni_id, uni_annot in unip.items():
            uni_name = uni_annot[self.UNIPROT_NAME]
            for mapping in uni_annot[self.MAPPINGS]:
                chain = mapping[self.CHAIN]
                if (pdb, chain) in self.pdbchain2uniidname and self.pdbchain2uniidname[(pdb, chain)] != (uni_id, uni_name):
                    old_id, old_name = self.pdbchain2uniidname[(pdb, chain)]
                    raise Exception(f'{pdb} chain {chain} has an ambiguous mapping to Uniprot: {old_id} vs {uni_id}')
                else:
                    self.pdbchain2uniidname[(pdb, chain)] = (uni_id, uni_name)


class RedirectIO:

    def __init__(self, stdin: Optional[str] = None, stdout: Optional[str] = None, stderr: Optional[str] = None, append_stdout: bool = False, append_stderr: bool = False):
        self.new_in_file = stdin
        self.new_out_file = stdout
        self.new_err_file = stderr
        self.append_stdout = append_stdout
        self.append_stderr = append_stderr

    def __enter__(self):
        if self.new_in_file is not None:
            self.new_in = open(self.new_in_file, 'r')
            self.old_in = sys.stdin
            sys.stdin = self.new_in
        if self.new_out_file is not None:
            self.new_out = open(self.new_out_file, 'a' if self.append_stdout else 'w')
            self.old_out = sys.stdout
            sys.stdout = self.new_out
        if self.new_err_file is not None:
            self.new_err = open(self.new_err_file, 'a' if self.append_stderr else 'w')
            self.old_err = sys.stderr
            sys.stderr = self.new_err

    def __exit__(self, exctype, excinst, exctb):
        if self.new_in_file is not None:
            sys.stdin = self.old_in
            self.new_in.close()
        if self.new_out_file is not None:
            sys.stdout = self.old_out
            self.new_out.close()
        if self.new_err_file is not None:
            sys.stderr = self.old_err
            self.new_err.close()


class ProgressBar:
    DONE_SYMBOL = '█'
    TODO_SYMBOL = '-'
    def __init__(self, n_steps, width=None, title='', prefix='', suffix='', writer=None):
        self.n_steps = n_steps # expected number of steps
        self.width = width if width is not None else shutil.get_terminal_size().columns
        self.width -= len(prefix) + len(suffix) + 12  # self.width -= len(prefix) + len(suffix) + 10
        self.title = (' '+title+' ')[0:min(len(title)+2, self.width)]
        self.prefix = prefix
        self.suffix = suffix
        self.writer = writer if writer is not None else sys.stdout
        self.done = 0 # number of completed steps
        self.shown = 0 # number of shown symbols

    def start(self):
        self.writer.write(' ' * len(self.prefix))
        self.writer.write(' ┌' + self.title + '─'*(self.width-len(self.title)) + '┐\n')
        self.writer.flush()
        self.step(0, force=True)
        return self

    def step(self, n_steps=1, force=False):
        if n_steps == 0 and not force:
            return
        self.done = min(self.done + n_steps, self.n_steps)
        progress = self.done / self.n_steps
        new_shown = int(self.width * progress)
        if new_shown != self.shown or force:
            self.writer.write(f'\r{self.prefix} └')
            self.writer.write(self.DONE_SYMBOL * new_shown + self.TODO_SYMBOL * (self.width - new_shown))
            self.writer.write(f'┘ {int(100*progress):>3}% {self.suffix} ')
            self.writer.flush()
            self.shown = new_shown  

    def finalize(self):
        self.step(self.n_steps - self.done)
        self.writer.write('\n')
        self.writer.flush()


Task = namedtuple('Task', ['name', 'function', 'args', 'kwargs', 'stdin', 'stdout', 'stderr', 'pre_message', 'post_message'])

class Pipeline:

    def __init__(self, checkpoint_file=None):
        self.checkpoint_file = checkpoint_file
        self.tasks = []
    
    def add_task(self, name: str, function: Callable, *args,
            stdin: Optional[str] = None, stdout: Optional[str] = None, stderr: Optional[str] = None, 
            pre_message: Optional[str] = None, post_message: Optional[str] = None, 
            kwargs = {}, **kwargs_):
        '''Append new task to the pipeline. 
        The task will be called as function(*args, **kwargs_). 
        If a keyword argument conflicts with a keyword arguments of add_task, use kwargs instead of kwargs_.
        Task names must be unique and not contain newline characters. 
        Task name can be None - such tasks are not logged into checkpoint file and will be repeated after resume.
        '''
        if name is not None and '\n' in name:
            raise ValueError(f'Task name must not contain newline character.')
        if name is not None and any(task.name == name for task in self.tasks):
            raise ValueError(f'Duplicate task name in pipeline: {name}. Use unique task names!')
        if len(kwargs) > 0 and len(kwargs_) > 0:
            raise ValueError('Use either kwargs or **kwargs_, not both!')
        kwargs_.update(kwargs)
        self.tasks.append(Task(name, function, args, kwargs_, stdin, stdout, stderr, pre_message, post_message))
    
    def start(self):
        '''Perform all tasks in the pipeline.'''
        if self.checkpoint_file is not None:
            with open(self.checkpoint_file, 'w', encoding=DEFAULT_ENCODING) as f:
                f.write('')
        for task in self.tasks:
            self._do_task(task)
    
    def resume(self):
        '''Perform all tasks in the pipeline, except for those completed earlier (i.e. listed in the checkpoint file).'''
        if self.checkpoint_file is not None:
            try: 
                with open(self.checkpoint_file, 'r', encoding=DEFAULT_ENCODING) as f:
                    completed_tasks = [line.rstrip('\n\r') for line in f]
            except OSError:
                completed_tasks = []
        else:
            completed_tasks = []
        for task in self.tasks:
            if task.name not in completed_tasks:
                self._do_task(task)
    
    def _do_task(self, task):
        if task.pre_message is not None:
            print(task.pre_message)
        with RedirectIO(stdin=task.stdin, stdout=task.stdout, stderr=task.stderr):
            result = task.function(*task.args, **task.kwargs)
        if self.checkpoint_file is not None and task.name is not None:
            with open(self.checkpoint_file, 'a', encoding=DEFAULT_ENCODING) as f:
                f.write(task.name + '\n')
        if task.post_message is not None:
            print(task.post_message)
