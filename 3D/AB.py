#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2023/01/21

import re
import sys
import math
import argparse
import numpy as np
import pandas as pd
from path import Path
from scipy import sparse
from collections import OrderedDict
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify, read_file, runshell


log = richlog()

def store_fpkm(fpkmdf):
    fpkmD = {}
    col_names = fpkmdf.columns.tolist()
    for i in range(1,fpkmdf.shape[1]):
        fpkmD[col_names[i]] = fpkmdf.set_index(col_names[0]).to_dict()[col_names[i]]
    return fpkmD


def assign_fpkm(beddf,fpkmD,prefix):
    for key in fpkmD.keys():
        N = 0
        concatL = []
        for index,row in beddf.iterrows():
            if not isinstance(row['ID'],float):
                fpkmL = []
                idL = row['ID'].split(',')
                for id in idL:
                    if id not in fpkmD[key].keys():
                        N +=1
                    else:
                        fpkmL.append(math.log2(fpkmD[key][id]+1))
                if len(fpkmL) == 0:
                    bin_mean = 0
                else:
                    bin_mean = sum(fpkmL)/len(fpkmL)
                new_row = pd.DataFrame({'seqname':[row['seqname']],
                'start':[row['start']],'end':[row['end']],
                'value':[bin_mean]})
                concatL.append(new_row)
            else:
                new_row = pd.DataFrame({'seqname':[row['seqname']],
                'start':[row['start']],'end':[row['end']],
                'value':np.nan})
                concatL.append(new_row)
        res = pd.concat(concatL)
        res.to_csv(f'{prefix}.{key}.FPKM.csv',header=False,index=False)
        log.info(f'Completed! the output file is {prefix}.{key}.FPKM.csv')
        print(f'Total missing number of sample {key} is {N}')


#outside command 
def CalBinFPKM(args):
    '''
    Calculate the average expression value of each bin
    >>> %(prog)s -f <fasta> -w windows_file [-o output] [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=CalBinFPKM.__name__,
                        description=CalBinFPKM.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('-e','--expmatrix',required=True,
            help='Input fasta file')
    pReq.add_argument('-b','--bed',required=True,
            help='Input the bin file containing the gene ID (available from utils.dist SlidingWindow2GeneID)')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.expmatrix)
    fpkmdf = read_file(args.expmatrix)
    beddf = read_file(args.bed,['seqname','start','end','ID'])
    prefix = re.findall('([\.\_0-9a-zA-Z]+)\.bed',args.bed)[0]
    fpkmD = store_fpkm(fpkmdf)
    assign_fpkm(beddf,fpkmD,prefix)
    

def main():
    actions = (
            ("CalBinFPKM", "Calculate the average expression value of each bin"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()