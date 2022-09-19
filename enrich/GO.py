#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/06/08


import re
import sys
import math
import argparse
import pandas as pd
from Bio import SeqIO
from path import Path
from natsort import natsorted
from collections import OrderedDict
from rich.traceback import install
from intervaltree import IntervalTree, Interval
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify, read_file
from PYHR.utils.fasta import Fa2dict
from PYHR.formats.gff import import_gff,get_gene_id

log = richlog()

def obo2db(file):
    num2row = OrderedDict()
    id2attr = OrderedDict()
    with open(file,'r',encoding='utf8') as f:
        for num,value in enumerate(f):
            num2row[num] = value.strip()
        for key in num2row.keys():
            if num2row[key].startswith('data-version'):
                ver = re.findall('data-version: (.*)',num2row[key])[0].split('/')[1]
            if '[Term]' in num2row[key]:
                goid = re.findall('id: (GO:[0-9]+)',num2row[key+1])[0]
                godes = re.findall('name: (.*)',num2row[key+2])[0]
                goclass = re.findall('namespace: (.*)',num2row[key+3])[0].replace('_',' ')
                # print(goid,godes,goclass,sep='\t',file=output)
                id2attr[goid] = [godes,goclass]
    with open('go.version_{}.txt'.format(ver),'w') as w:
        w.write('GO\tDescription\tlevel\n')
        for i in id2attr.keys():
            w.write('{}\t{}\t{}\n'.format(i,id2attr[i][0],id2attr[i][1]))
    return ver


#outside command 
def convertOBO(args):
    '''
    Convert go-basic.obo to an easy-to-read file
    >>> %(prog)s <go-basic.obo> [output] [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=convertOBO.__name__,
                        description=convertOBO.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')


    pReq.add_argument('obo',
            help='Input go-basic obo file')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.obo)
    ver = obo2db(args.obo)
    log.info('Completed! the output file is `go.version_{}.txt`'.format(ver))


def main():
    actions = (
            ("convertOBO", "Convert go-basic.obo to an easy-to-read file"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
