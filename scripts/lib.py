import sys
import json
import requests
from typing import List, Dict, Tuple, Union, Iterator, Iterable, Callable, TypeVar

from constants import *

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
        with open(conversion_table_file) as r:
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
