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

class SkipException(Exception):
    pass

class Skippo:
    def __enter__(self):
        raise SkipException('Unicorn')
    def __exit__(self, extype, exval, extrace):
        # if extype == SkipException:
        print(f'Skipped because of {exval}')
        return True

def pica(first, *args, stdout=None, **kwargs):
    print('FIRST:', first)
    print('ARGS:', args)
    print('KWARGS', kwargs)

#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
    print(vars(args))
    return vars(args)


def main() -> Optional[int]:
    pica('januar', 'teply', 'vlhky', temp='15C', stdout='out')
    print()
    pica('januar', 'teply', 'vlhky', kwargs={'temp':'15C'})
    x = 5
    f = lambda: print(x)
    x = 10
    f()
    # with Skippo():
    #     print('Never')


if __name__ == '__main__':
    args = parse_args()
    exit_code = main(**args)
    if exit_code is not None:
        exit(exit_code)