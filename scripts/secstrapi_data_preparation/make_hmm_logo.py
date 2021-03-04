'''
This Python3 script does foo ...

Example usage:
    python3  foo.py  --foo 4  foo.txt 
'''
# TODO add description and example usage in docstring

import argparse
import numpy as np
import pandas as pd
from typing import Dict, Any, Optional, Tuple, Union

try:
    import no_gap_align
except ImportError:
    from . import no_gap_align

#  CONSTANTS  ################################################################################

AMINOACIDS = list('ACDEFGHIKLMNPQRSTVWY')

#  FUNCTIONS  ################################################################################


#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('hmmlogo_file', help='Output of running "hmmlogo *.hmm"', type=str)
    parser.add_argument('output_file', help='Filename for output image (PNG/TIF)"', type=str)
    parser.add_argument('--positions', help=f'Residue range from:to (1-based, to excluded, default = ":" (meaning start:end))', type=str, default=':')
    parser.add_argument('--title', help=f'Title for the logo', type=str, default='')
    parser.add_argument('--shift_numbers', help=f'Shift residue numbers by this number (default: 0)', type=int, default=0)
    # TODO add command line arguments
    args = parser.parse_args()
    return vars(args)


def main(hmmlogo_file: str, output_file: str, positions: Union[Tuple[Optional[int], Optional[int]], str] = ':', title: str = '', shift_numbers: Union[int, str] = 0) -> Optional[int]:
    # TODO add parameters
    '''Foo'''
    # TODO add docstring
    if isinstance(positions, str):
        positions = parse_range(positions)
    fro1, to1 = positions
    fro0, to0 = range_1_to_0_based(positions)
    shift_numbers = int(shift_numbers)
        
    heights, occupancy = parse_hmmlogo_file(hmmlogo_file)
    heights = heights[fro0:to0, :]
    # print(heights / heights.sum(axis=1, keepdims=True))
    occupancy = occupancy[fro0:to0]
    heights = pd.DataFrame(heights, columns=AMINOACIDS)
    first_index = fro1 if fro1 is not None else 1
    print(heights)
    print(occupancy)
    no_gap_align.run_logomaker_from_matrix(heights, occupancy, output_file, first_index=first_index+shift_numbers, title=title)

def parse_range(range_string: str) -> Tuple[int, int]:
    sfro, sto = range_string.split(':')
    fro = int(sfro) if sfro.strip() != '' else None
    to = int(sto) if sto.strip() != '' else None
    return (fro, to)

def range_1_to_0_based(the_range: Tuple[Optional[int], Optional[int]]) -> Tuple[Optional[int], Optional[int]]:
    fro, to = the_range
    if fro is not None:
        fro -=1
    if to is not None:
        to -=1
    return (fro, to)

def parse_hmmlogo_file(filename: str):
    with open(filename) as r:
        lines = [line.strip() for line in r.readlines()]
    lines = [line for line in lines if line != '']
    start_res_heights = lines.index('Residue heights')
    start_indels = lines.index('Indel values')
    res_height_lines = lines[start_res_heights+1:start_indels]
    indel_lines = lines[start_indels+1:]
    res_heights = []
    for line in res_height_lines:
        str_values = line.split(':')[1].split('(')[0].split()
        values = [float(x) for x in str_values]
        res_heights.append(values)
    indels = []
    for line in indel_lines:
        str_values = line.split(':')[1].split()
        values = [float(x) for x in str_values]
        indels.append(values)
    n = len(res_heights)
    res_heights = np.array(res_heights)
    indels = np.array(indels)
    assert res_heights.shape == (n, 20)
    assert indels.shape == (n, 3)
    occupancy = indels[:, -1]
    return res_heights, occupancy
            


if __name__ == '__main__':
    args = parse_args()
    exit_code = main(**args)
    if exit_code is not None:
        exit(exit_code)