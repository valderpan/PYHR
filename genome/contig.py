#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/05/10


import argparse
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.utils.fasta import Fa2dict


log = richlog()

def cal_ctgN50(CtgD,Nratio):
    Nxratio = [int(i) for i in Nratio[0].split(',')]
    ctgs_length = 0
    ctgs_num = 0
    for key in CtgD.keys():
        ctgs_num += 1
        ctgs_length += len(CtgD[key])

    log.info('Total number of ctgs:  {}'.format(ctgs_num))
    log.info('Total length of ctgs:  {} bp'.format(ctgs_length))
    log.info('Average contig size:  {} bp'.format(round(ctgs_length / ctgs_num, 2)))

    log.info('Longest ctg sequence length:  {} bp'.format(
        len(list(enumerate(sorted(CtgD.items(), key=lambda x: len(x[1]))))[-1][1][1])))
    log.info('Shortest ctg sequence length:  {} bp'.format(
        len(list(enumerate(sorted(CtgD.items(), key=lambda x: len(x[1]))))[0][1][1])))

    for nx in Nxratio:
        if nx > 100 or nx < 0:
            log.error('The nx value must be set bwteen 0 to 100 !')
        else:
            NLen = ctgs_length * (nx/100)
            N = 0
            for num,key in enumerate(sorted(CtgD.items(), key=lambda x: len(x[1]), reverse=True)):
                N += len(CtgD[key[0]])
                if N >= NLen:
                    log.info('N{} contig:  {} ({} sequences),length  {} bp'.format(nx,key[0],num+1 ,len(CtgD[key[0]])))
                    break


#outside command 
def N50(args):
    '''
    Calculate N50 for fsata
    >>> %(prog)s <fasta> <-nx nxratio> [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=N50.__name__,
                        description=N50.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fasta', 
            help='Input fasta file')
    pReq.add_argument('-nx','--nxratio', nargs='+', 
            help="Nx values to be printed seperated by ',' e.g. 50 for N50, 75 for N75")
    pOpt.add_argument('-h', '--help', action='help',
            help='Show help message and exit.')
    
    args = p.parse_args(args)

    check_file_exists(args.fasta)
    ctgD = Fa2dict(args.fasta)
    cal_ctgN50(ctgD,args.nxratio)


def main():
    actions = (
            ("N50", "Calculate N50(or N60,N70,...) for fsata"),
            # ("ExtractFasta", "Extract query seq by header"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()