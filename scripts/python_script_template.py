'''
This Python3 script does foo ...

Example usage:
    python3  foo.py  --foo 4  foo.txt 
'''
# TODO add description and example usage in docstring

import argparse
from typing import Dict, Any, Optional

#  CONSTANTS  ################################################################################


#  FUNCTIONS  ################################################################################


#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('foo', help='Dummy argument foo', type=str)
    parser.add_argument('--bar', help=f'Dummy option bar (default = {123})', type=int, default=123)
    parser.add_argument('--baz', help=f'Dummy switch baz', action='store_true')
    # TODO add command line arguments
    args = parser.parse_args()
    return vars(args)


def main(foo: str, bar: int = 123, baz: bool = False) -> Optional[int]:
    # TODO add parameters
    '''Foo'''
    # TODO add docstring
    pass
    # TODO add implementation


if __name__ == '__main__':
    args = parse_args()
    exit_code = main(**args)
    if exit_code is not None:
        exit(exit_code)