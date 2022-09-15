#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/05/08


import re
import sys
import argparse
from Bio import SeqIO
from rich import print
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify


log = richlog()

def Fa2dict(fasta_file):
    fastaD = {}
    for seq in SeqIO.parse(fasta_file,'fasta'):
        fastaD[seq.id] = seq.seq
    return fastaD


class FastaFilter():
    def __check_keywords(func):
        def inner(self,keywords):
            if not isinstance(keywords,str):
                raise TypeError('Please pass a string for keywords')
            return func(self,keywords)
        return inner
    
    def __init__(self,fasta):
        self.fastaD = Fa2dict(fasta)

    @__check_keywords
    def filter_start(self,keywords):
        newD = {key:self.fastaD[key] for key in self.fastaD if key.startswith(keywords)}
        return newD
    
    @__check_keywords
    def filter_end(self,keywords):
        newD = {key:self.fastaD[key] for key in self.fastaD if key.endswith(keywords)}
        return newD
    
    @__check_keywords
    def filter_re(self,keywords):
        newD = {key:self.fastaD[key] for key in self.fastaD if re.findall(keywords, key)[0] in key}
        return newD
    
    @__check_keywords
    def filter_match(self,keywords):
        newD = {key:self.fastaD[key] for key in self.fastaD if key == keywords}
        return newD
    
    @staticmethod
    def outputfastaD(fastaD,output):
        for key in fastaD.keys():
            print(key,file=output)
            print(fastaD[key],file=output)


def read_seqID(file):
    with open(file,'r') as f:
        seqIDL = [l for l in (line.strip() for line in f)]
    log.info('{} seqs are obtained'.format(len(seqIDL)))
    return seqIDL


def ExtractandOut_seq(fastaD,ID,output):
    unmatch = []
    match = 0
    with open(output,'w') as w:
        for key in fastaD.keys():
            if key in ID:
                match += 1
                w.write(">{}\n".format(key))
                w.write(str(fastaD[key])+'\n')
            else:
                unmatch.append(key)
    log.info('Total {} seqs are matched'.format(match))
    if len(unmatch)>0:
        if len(unmatch) > 3:
            log.info('Total {} seqs did are not matched, such as {}, {}, {},...'.format(len(unmatch),unmatch[0],unmatch[1],unmatch[2]))
        else:
            log.info('Only {} seqs did are not matched, such as {}...'.format(len(unmatch),unmatch[0]))


def CalEGS(fastaD):
    egsL = [];total_len = [];total_A = [];total_C = [];total_G = [];total_T = [];total_N = []
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('#seq','len','A','C','G','T','N'))
    for key in fastaD.keys():
        seq = str(fastaD[key])
        seqlen = len(seq)
        A = seq.upper().count('A')
        T = seq.upper().count('T')
        G = seq.upper().count('G')
        C = seq.upper().count('C')
        N = seq.upper().count('N')
        print(key,seqlen,A,C,G,T,N,sep='\t')
        egsL.extend([A,C,G,T])
        total_A.extend([A])
        total_C.extend([C])
        total_C.extend([G])
        total_C.extend([T])
        total_N.extend([N])
        total_len.extend([seqlen])
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('Total',sum(total_len),sum(total_A),
                                        sum(total_C),sum(total_G),sum(total_T),sum(total_N)))
    print('Effective genome size: {}'.format(sum(total_len)-sum(total_N)))

def main():
    actions = (
            ("FilterFasta", "Filter fasta by header"),
            ("ExtractFasta", "Extract query seq by header"),
            ('CalEffGenomeSize', "Calculate the effective size of the genome")
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


#outside command 
def FilterFasta(args):
    '''
    Filter fasta by header
    >>> %(prog)s <fasta> <-m mode> <-k keywords> <-o output> [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=FilterFasta.__name__,
                        description=FilterFasta.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fasta', 
            help='Input fasta file')
    pReq.add_argument('-m','--mode', choices=['start','end','re','match'], 
            help='Choose filrter mode')
    pReq.add_argument('-k','--keywords', 
            help='Input the filtered keywords or Regular expressions')
    pReq.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='Output file [default: stdout]')
    
    pOpt.add_argument('-h', '--help', action='help',
            help='Show help message and exit.')
    
    args = p.parse_args(args)

    check_file_exists(args.fasta)
    f = FastaFilter(args.fasta)
    if not isinstance(args.mode,str):
        raise TypeError('Please pass a string for mode')
    else:
        if args.mode not in  ['start','end','re','match']:
            raise ValueError('The input mode must be one of the following four :[start,end,re,match]')

    if args.mode == 'start':        
        newD = f.filter_start(args.keywords)
    elif args.mode == 'end':
        newD = f.filter_end(args.keywords)
    elif args.mode == 'match':
        newD = f.filter_match(args.keywords)
    else:
        newD = f.filter_re(args.keywords)

    f.outputfastaD(newD,args.output)


def ExtractFasta(args):
    '''
    Filter fasta by header
    >>> %(prog)s <fasta> <-q query_seqID file> <-o output> [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=ExtractFasta.__name__,
                        description=ExtractFasta.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fasta', 
            help='Input fasta file')
    pReq.add_argument('-q','--query', 
            help='Input query seqid file')
    pReq.add_argument('-o', '--output',
            help='Output file')
    pOpt.add_argument('-h', '--help', action='help',
            help='Show help message and exit.')
    
    args = p.parse_args(args)

    check_file_exists(args.fasta)
    fastaD = Fa2dict(args.fasta)
    check_file_exists(args.query)
    seqIDL = read_seqID(args.query)
    
    ExtractandOut_seq(fastaD,seqIDL,args.output)


def CalEffGenomeSize(args):
    '''
    Calculate the effective size of the genome
    >>> %(prog)s <fasta> <-o output> [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=ExtractFasta.__name__,
                        description=ExtractFasta.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fasta', 
            help='Input fasta file')
    pReq.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='Output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='Show help message and exit.')
    
    args = p.parse_args(args)
    check_file_exists(args.fasta)
    fastaD = Fa2dict(args.fasta)
    CalEGS(fastaD)


if __name__ == "__main__":
    main()