#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/05/18

import re
import sys
import argparse
from rich import print
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify


log = richlog()

def read_result(file):
    rowL = []
    with open(file) as f:
        lines = (line.strip() for line in f)
        for line in lines:
            rowL.append(line)
    return rowL


def re_findends(patten,rowL):
    for i in rowL:
        if i.endswith(patten):
            return int(re.findall('([0-9]+)',i)[0])
        else:
            continue


def parse_result(rowL):
    single_reads = re_findends('reads; of these:',rowL)
    uniq_reads1 = re_findends('aligned concordantly exactly 1 time',rowL)
    uniq_reads2 = re_findends('aligned discordantly 1 time',rowL)
    uniq_reads3 = re_findends('aligned exactly 1 time',rowL)
    multi_reads1 = re_findends('aligned concordantly >1 times',rowL)
    multi_reads2 = re_findends('aligned >1 times',rowL)
    return [single_reads,uniq_reads1,uniq_reads2,uniq_reads3,multi_reads1,multi_reads2]


def cal_rate(reads_num_list):
    Total_reads = reads_num_list[0]*2
    uniq_reads = reads_num_list[1]*2+reads_num_list[2]*2+reads_num_list[3]
    uniq_rate = uniq_reads/Total_reads*100
    multi_reads = reads_num_list[4]*2+reads_num_list[5]
    multi_rate = multi_reads/Total_reads*100
    return [Total_reads,uniq_reads,uniq_rate,multi_reads,multi_rate]


def print_output(logfile,resL,output=None):
    print('Sample file : {}'.format(logfile),file=output)
    print('    Total reads num : {}'.format(resL[0]),file=output)
    print('    Uniq reads num : {}'.format(resL[1]),file=output)
    print('    Uniq reads rate : {}%'.format(round(resL[2],2)),file=output)
    print('    Multi reads num : {}'.format(resL[3]),file=output)
    print('    Multi reads rate : {}%'.format(round(resL[4],2)),file=output)



#outside command 
def StatHisat2MappedRate(args):
    '''
    Calculate the comparison of each indicator after hisat2 runs by specifying the results obtained with the --summary-file parameter
    >>> %(prog)s <stat.file> [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=StatHisat2MappedRate.__name__,
                        description=StatHisat2MappedRate.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')


    pReq.add_argument('stat',
            help='Input the statistics file specified by hisat2 --summary-file')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            #default=sys.stdout, 
            help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.stat)
    rowL = read_result(args.stat)
    readsL = parse_result(rowL)
    resL = cal_rate(readsL)

    if args.output :
        print_output(args.stat,resL,args.output)
        log.info('Completed! The output file is `{}` '.format(args.output.name))
    else:
        print_output(args.stat,resL)


def main():
    actions = (
            ("StatHisat2MappedRate", "Calculate the comparison of each indicator after hisat2 runs by specifying the results obtained with the --summary-file parameter"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
