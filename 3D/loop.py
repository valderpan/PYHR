#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/10/08


import re
import sys
import argparse
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

def read_loop(file):
    '''
    The loop is generated by Hicup + Hichipper +filter_siginificant_loop pipeline
    The basic format of a loop file is as follows:

    X1	X2	X3	X4	X5	X6	X7	X8	size	corrected_strength	loop_fdr
    chr10	102800185	102803154	chr10	102825482	102828931	.	3	6418	2	9.012691140774453e-5
    chr10	120957506	120959815	chr10	121204263	121206068	.	3	4114	3	4.177403991386223e-9
    '''
    loopdf = read_file(file)
    loopdf = loopdf.sort_values(['X1','X2'],ignore_index=True)
    return loopdf


def loopintervale(loopdf):
    loopinter = {}
    for index,row in loopdf.iterrows():
        if not row['X1'] in loopinter:
            loopinter[row['X1']] = IntervalTree()
            loopinter[row['X1']].addi(int(row['X2']),int(row['X3']),
            {'{}-{}'.format(row['X1'],index):'{}-{}-{}-{}'.format(row['X4'],row['X5'],row['X6'],index)})
        else:
            loopinter[row['X1']].addi(int(row['X2']),int(row['X3']),
            {'{}-{}'.format(row['X1'],index):'{}-{}-{}-{}'.format(row['X4'],row['X5'],row['X6'],index)})
    return loopinter


def parse_inter(loop_i1,loop_i2):
    loop_i1_q = [];loop_i2_q =[]
    for chr in loop_i1.keys():
        if chr in loop_i2.keys():
            for inter in sorted(loop_i1[chr]):
                if loop_i2[chr].overlap(inter):
                    other_a1 = list(inter.data.values())[0]
                    chr1,start1,end1,index1 = other_a1.split('-')
                    a1_inter = Interval(int(start1),int(end1),'{}-{}'.format(chr1,index1))
                    for i in loop_i2[chr].overlap(inter):
                        other_a2 = list(i.data.values())[0]
                        chr2,start2,end2,index2 = other_a2.split('-')
                        a2_inter = Interval(int(start2),int(end2),'{}-{}'.format(chr1,index1))
                        if a1_inter.overlaps(a2_inter):
                            loop_i1_q.append(int(index1))
                            loop_i2_q.append(int(index2))
    return loop_i1_q,loop_i2_q

def get_querydf(index,refdf):
    qdf = refdf[refdf.index.isin(index)]
    qdf = qdf.sort_values(['X1','X2'])
    return qdf


def get_querydf_other(index,refdf):
    qdf = refdf[~refdf.index.isin(index)]
    qdf = qdf.sort_values(['X1','X2'])
    return qdf


## outside command 
def findOverlapLoops(args):
    """
    Find loops with overlap in both left and right anchors in both sets of data
    >>> %(prog)s -a loop_a -b loop_b [Options]
    Output:
    .overlap.bedpe: loops having overlap with another set of data
    .specific.bedpe: loops that no overlap with another set of data, and it is self-specific loop
    """
    install()
    p = argparse.ArgumentParser(prog=findOverlapLoops.__name__,
                        description=findOverlapLoops.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-a','--loopa', 
            help='Input the first set of loop data')
    pReq.add_argument('-b', '--loopb',default=None,
            help='Input the second set of loop data')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    check_file_exists(args.loopa)
    check_file_exists(args.loopb)
    loopadf = read_loop(args.loopa)
    loopbdf = read_loop(args.loopb)
    log.info('Save {} loop data to IntervalTree...'.format(args.loopa))
    loopa_inter = loopintervale(loopadf)
    log.info('Save {} loop data to IntervalTree...'.format(args.loopb))
    loopb_inter = loopintervale(loopbdf)
    log.info('Find overlap Loops....')
    loopa_q,loopb_q = parse_inter(loopa_inter,loopb_inter)

    
    common_adf = get_querydf(loopa_q,loopadf)
    specific_adf = get_querydf_other(loopa_q,loopadf)
    log.info('The number of original loops of {} is {}'.format(args.loopa,loopadf.shape[0]))
    log.info('The number of common loops of {}.overlap.bedpe is {}'.format(args.loopa,common_adf.shape[0]))
    log.info('The percentage of overlap loops is {}'.format(round(common_adf.shape[0]/loopadf.shape[0],2)*100))
    log.info('The number of specific loop of {} is {}'.format(args.loopa,loopadf.shape[0]-common_adf.shape[0]))
    common_adf.to_csv('{}.overlap.bedpe'.format(args.loopa),sep='\t',header=False,index=False)
    specific_adf.to_csv('{}.specific.bedpe'.format(args.loopa),sep='\t',header=False,index=False)


    common_bdf = get_querydf(loopb_q,loopbdf)
    specific_bdf = get_querydf_other(loopb_q,loopbdf)
    log.info('The number of original loops of {} is {}'.format(args.loopb,loopbdf.shape[0]))
    log.info('The number of common loops of {}.overlap.bedpe is {}'.format(args.loopb,common_bdf.shape[0]))
    log.info('The percentage of overlap loops is {}'.format(round(common_bdf.shape[0]/loopbdf.shape[0],2)*100))
    log.info('The number of specific loop of {} is {}'.format(args.loopb,loopbdf.shape[0]-common_bdf.shape[0]))
    common_bdf.to_csv('{}.overlap.bedpe'.format(args.loopb),sep='\t',header=False,index=False)
    specific_bdf.to_csv('{}.specific.bedpe'.format(args.loopb),sep='\t',header=False,index=False)

    
def main():
    actions = (
            ("findOverlapLoops", "Find loops with overlap in both left and right anchors in both sets of data"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()