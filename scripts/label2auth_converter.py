import sys
from typing import Tuple

class Label2AuthConverter:
    def __init__(self, conversion_table_file: str):
        self.file = conversion_table_file
        self.table = {}
        self.chain_table = {}
        with open(conversion_table_file) as r:
            for line in r:
                if not line.lstrip().startswith('#'):
                    chain, resi, auth_chain, auth_resi, auth_ins_code = line.strip().split('\t')
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
            raise Exception(f'Did not find residue {chain} {resi} in {self.file}')

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
