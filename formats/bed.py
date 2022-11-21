#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/05/08


"""
bed format and some util tools.
"""


import sys
import argparse
import pandas as pd
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify


log = richlog()

def main():
    actions = (
            ("RandomBed", "random select several gene from bed file."),
            ("CountBedBySeq", "count bed file by seqid"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def import_bed6(bed, sort=False):
    """
    import 6-columns bed file
    """
    names = ['seqid', 'start', 'end', 'name', 'score', 'strand']
    check_file_exists(bed)
    df = pd.read_csv(bed, sep='\t', header=None, index_col=None, names=names)
    if sort:
        df = df.sort_values(by=['seqid', 'start'])
    return df


## out command
def CountBedBySeq(args):
    """
    count bed file by seqid
    >>> %(prog)s <bed> [Options]
    """
    install()
    p = argparse.ArgumentParser(prog=CountBedBySeq.__name__,
                        description=CountBedBySeq.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bed', 
            help='input bed file')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    df = import_bed6(args.bed)
    df = df.groupby('seqid').size().reset_index(name='counts')
    df.to_csv(args.output, sep='\t', line_terminator="\n", header=False, index=False)


# TODO 下面这个函数可以有但没必要
def RandomBed(args):
    """
    random select several genes from bed file.
    >>> %(prog)s <gene.bed> <--target target.bed/-n nums> [Options]
    """
    install()
    p = argparse.ArgumentParser(prog=__file__,
                        description=RandomBed.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bed', help='source bed file')
    pOpt.add_argument('--target', help='target bed file')
    pOpt.add_argument('-n', '--nums', type=int, 
            help='select numbers (will be overlaped with `--target`)')
    pOpt.add_argument('-o', '--out', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('--include', action='store_true', default=False,
            help='if include the target gene in source bedfile [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    source_df = pd.read_csv(args.bed, sep='\t', header=None,
            names=('seqid', 'start', 'end', 'gene'))
    if not args.target:
        if args.nums:
            nums = args.nums
        else:
            log.error('Must input `--nums` or `--target` options')
    else:
        target_df = pd.read_csv(args.target, sep='\t', header=None,
        names=('seqid', 'start', 'end', 'gene'))
        nums = len(target_df)

    if args.target and not args.include:
        source_df = source_df[~source_df['gene'].isin(target_df['gene'])]

    random_df = source_df.sample(nums)
    random_df = random_df.sort_values(by=['seqid', 'start'])
    random_df.to_csv(args.out, sep='\t', header=False,index=False)


if __name__ == "__main__":
    main()