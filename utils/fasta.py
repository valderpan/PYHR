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
from rich.console import Console
from rich.table import Column, Table
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify,read_file


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
        with open(output,'w') as w1:
            for key in fastaD.keys():
                w1.write('>'+key+'\n')      
                w1.write(str(fastaD[key])+'\n')


def read_seqID(file):
    with open(file,'r') as f:
        seqIDL = [l for l in (line.strip() for line in f)]
    log.info('{} seqs are obtained'.format(len(seqIDL)))
    return seqIDL


def ExtractandOut_seq(fastaD,ID,output):
    unmatch = [];obm = [];matched = []
    with open(output,'w') as w:
        for i in ID:
            if i in fastaD.keys():
                matched.append(i)
                w.write(">{}\n".format(i))
                w.write(str(fastaD[i])+'\n')
            else:
                unmatch.append(i)
    if len(unmatch) > 0:
        if len(unmatch) > 3:
            log.warning('Total {} seqs did are not matched, such as {}, {}, {},...'.format(len(unmatch),unmatch[0],unmatch[1],unmatch[2]))
        else:
            log.warning('Only {} seqs did are not matched, such as {}...'.format(len(unmatch),unmatch))


def CalEGS(fastaD):
    console = Console()
    table = Table(show_header=True, header_style="bold magenta")
    egsL = [];total_len = [];total_A = [];total_C = [];total_G = [];total_T = [];total_N = []
    table.add_column("Seq",width=12)
    table.add_column("Len")
    table.add_column("A")
    table.add_column("T")
    table.add_column("G")
    table.add_column("C")
    table.add_column("N")
    # print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('#seq','len','A','C','G','T','N'))
    for key in fastaD.keys():
        seq = str(fastaD[key])
        seqlen = len(seq)
        A = seq.upper().count('A')
        T = seq.upper().count('T')
        G = seq.upper().count('G')
        C = seq.upper().count('C')
        N = seq.upper().count('N')
        # print(key,seqlen,A,C,G,T,N,sep='\t')
        table.add_row(key, str(seqlen), str(A), str(T), str(G), str(C), str(N))
        egsL.extend([A,C,G,T])
        total_A.extend([A])
        total_C.extend([C])
        total_G.extend([G])
        total_T.extend([T])
        total_N.extend([N])
        total_len.extend([seqlen])
    # print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('Total',sum(total_len),sum(total_A),
    #                                     sum(total_C),sum(total_G),sum(total_T),sum(total_N)))
    table.add_row('Total', str(sum(total_len)), 
        str(sum(total_A)), str(sum(total_T)), str(sum(total_G)), str(sum(total_C)), str(sum(total_N)))
    table.add_row('[red]Effective genome size[/red]', '[red]{}[/red]'.format(str(sum(total_len)-sum(total_N))))
    # print('Effective genome size: {}'.format(sum(total_len)-sum(total_N)))
    console.print(table)


def replace_fastaheader(fastaD,file,output):
    df = read_file(file,['oldid','newid'])
    old2new = df.set_index('oldid').to_dict()['newid']
    if len(old2new) < len(fastaD):
        log.warning('The content provided in the reference file does not contain all the headers,those that are not included will be printed as ctg_*')
    with open(output,'w') as w:
        for key in fastaD.keys():
            if key in old2new:
                w.write('>{}\n'.format(old2new[key]))
                w.write('{}\n'.format(fastaD[key]))
            else:
                w.write('>ctg_{}\n'.format(key))
                w.write('{}\n'.format(fastaD[key]))


def main():
    actions = (
            ("FilterFasta", "Filter fasta by header"),
            ("ExtractFasta", "Extract query seq by header"),
            ('CalEffGenomeSize', "Calculate the effective size of the genome"),
            ('ReplaceFastaHeader',"Replace the original header in the fasta file")
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
    pReq.add_argument('-o', '--output',
            help='Output file')
    pOpt.add_argument('-h', '--help', action='help',
            help='Show help message and exit.')
    
    args = p.parse_args(args)

    check_file_exists(args.fasta)
    f = FastaFilter(args.fasta)
    if not isinstance(args.mode,str):
        raise TypeError('Please pass a string for mode')
    else:
        if not args.mode in  ['start','end','re','match']:
            raise ValueError('The input mode must be one of the following four :[start,end,re,match]')

    if args.mode == 'start':        
        newD = f.filter_start(args.keywords)
    elif args.mode == 'end':
        newD = f.filter_end(args.keywords)
    elif args.mode == 'match':
        newD = f.filter_match(args.keywords)
    else:
        newD = f.filter_re(args.keywords)
    if len(newD) >=1 :
        log.info('A total of {} seqs were saved after filtering'.format(len(newD)))
        f.outputfastaD(newD,args.output)
        log.info('The results file is output to `{}` '.format(args.output))
    else:
        log.error('None of the seqs are saved after filtering, please check the matching conditions! ')
    

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
    # pReq.add_argument('-o', '--output', type=argparse.FileType('w'),
    #         default=sys.stdout, help='Output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='Show help message and exit.')
    
    args = p.parse_args(args)
    check_file_exists(args.fasta)
    fastaD = Fa2dict(args.fasta)
    CalEGS(fastaD)


def ReplaceFastaHeader(args):
    '''
    Replace the original header in the fasta file
    ==============================================================================================================
    Note:The input replacement reference file must be in two-column format and have no headers !
        eg: GK000010.2   chr10
            GK000011.2   chr11
    Make sure that all pseudochromosome ids are in the reference file, as unincluded ids will be output as ctg_* !
    ==============================================================================================================
    >>> %(prog)s -i <fasta> -r <replace ref> -o <output fasta> [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=ReplaceFastaHeader.__name__,
                        description=ReplaceFastaHeader.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fasta', 
            help='Input fasta file')
    pReq.add_argument('reference', 
            help='Input replance reference file, it must be in two-column format and have no headers')
    pReq.add_argument('output',
            help='Output fasta file')
    pOpt.add_argument('-h', '--help', action='help',
            help='Show help message and exit.')
    
    args = p.parse_args(args)
    check_file_exists(args.fasta)
    fastaD = Fa2dict(args.fasta)
    check_file_exists(args.reference)
    replace_fastaheader(fastaD,args.reference,args.output)
    log.info('Completed! the output file is `{}`'.format(args.output))

if __name__ == "__main__":
    main()