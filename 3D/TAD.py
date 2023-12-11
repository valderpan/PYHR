#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2023/12/04


import re
import sys
import argparse
import numpy as np
import pandas as pd
from path import Path
from natsort import natsorted
from collections import OrderedDict
from rich.traceback import install
from intervaltree import IntervalTree,Interval
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify, read_file
from PYHR.utils.fasta import Fa2dict

log = richlog()


def parseBoundaries(Boundaries_bed):
    with open(Boundaries_bed,'r') as f:
        prefix = re.findall('([0-9A-Za-z\.\_]+)\-\-',Boundaries_bed)[0]
        with open(f'{prefix}.TAD.bed','w') as w:
            lines = f.readlines()
            for i in range(0,len(lines)):
                if i == 0 or i+1 >= len(lines):
                    continue
                else:
                    fields1 = lines[i].strip().split('\t')
                    fields2 = lines[i+1].strip().split('\t')
                    w.write(f"{fields1[0]}\t{fields1[2]}\t{fields2[1]}\n")
    return prefix


## outside command
def boundaries2bed(args):
    """
    Convert boundaries to bed based on cworld results
    >>> %(prog)s <bed> [iced.matrix][Options]
    """ 
    install()
    p = argparse.ArgumentParser(prog=boundaries2bed.__name__,
                        description=boundaries2bed.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    
    pReq.add_argument('BoundariesBed',
                    type=str,
                    help="Input cworld's result file boundaries bed")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.BoundariesBed)
    res_prefix = parseBoundaries(args.BoundariesBed)
    log.info(f'The results file is output to `{res_prefix}.TAD.bed`')


def main():
    actions = (
            ("boundaries2bed", "Convert boundaries to bed based on cworld results"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()