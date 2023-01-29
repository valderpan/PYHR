#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/06/02


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
from PYHR.formats.gff import import_gff,get_Attr_dict,get_gene_id,get_gene_name


log = richlog()

def stat_seqLen(fastaD):
    fasta2Len = {}
    for key in fastaD.keys():
        fasta2Len[key] = len(fastaD[key])
        print('{}\t{}'.format(key,len(fastaD[key])))
    return fasta2Len


def sliding_windows(fastaD,window,step,output):
    with open(output,'w') as w:
        for key in fastaD.keys():
            start = 1
            while start < fastaD[key]:
                if start+window < fastaD[key]:
                    w.write('{}\t{}\t{}\n'.format(key,start,start+window-1))
                    start += step
                else:
                    w.write('{}\t{}\t{}\n'.format(key,start,fastaD[key]))
                    break


def gff2tree_id(gff3,type):
    chr2tree = {}
    gff = import_gff(gff3,'gene')
    geneid = []
    for _,row in gff.iterrows():
        db = get_Attr_dict(row['attributes'])
        if type == 'ID':
            geneid.append(get_gene_id(db))
        elif type == 'Name':
            geneid.append(get_gene_name(db))
    gff['id'] = geneid
    gff = gff.loc[:,['start','end','id']].reset_index()
    for _,row in gff.iterrows():
        if row['seqid'] not in chr2tree:
            chr2tree[row['seqid']] = IntervalTree()
            chr2tree[row['seqid']].addi(int(row['start']),int(row['end']),row['id'])
        else:
            chr2tree[row['seqid']].addi(int(row['start']),int(row['end']),row['id'])
    return chr2tree


def gff2tree_len(gff3):
    chr2tree = {}
    gff = import_gff(gff3,'gene')
    geneid = []
    for _,row in gff.iterrows():
        db = get_Attr_dict(row['attributes'])
        geneid.append(get_gene_id(db))
    gff['id'] = geneid
    gff = gff.loc[:,['start','end','id']].reset_index()
    for _,row in gff.iterrows():
        if row['seqid'] not in chr2tree:
            chr2tree[row['seqid']] = IntervalTree()
            chr2tree[row['seqid']].addi(int(row['start']),int(row['end']),int(row['end']-int(row['start'])+1))
        else:
            chr2tree[row['seqid']].addi(int(row['start']),int(row['end']),int(row['end']-int(row['start'])+1))
    return chr2tree


def window2tree(windows_file):
    win2tree = {}
    with open(windows_file) as f:
        lines = (line.strip() for line in f)
        for line in lines:
            line_list = line.split()
            if line_list[0] not in win2tree:
                win2tree[line_list[0]] = IntervalTree()
                win2tree[line_list[0]].addi(int(line_list[1]),int(line_list[2]))
            else:
                win2tree[line_list[0]].addi(int(line_list[1]),int(line_list[2]))
    return win2tree


def statGeneNum(chr_tree,wintree):
    for key in wintree.keys():
        if 'ctg' not in key and 'tig' not in key and 'utg' not in key:
            for inter in sorted(wintree[key]):
                all_cont = 0
                if chr_tree[key].overlap(inter):

                    for i in chr_tree[key].overlap(inter):
                        if i.begin >= inter.begin and i.end < inter.end:
                            all_cont +=1
                        elif i.end >= inter.begin and i.begin < inter.begin:
                            if i.end-inter.begin >= (i.end-i.begin)/2:
                                all_cont += 1

                            else:
                                continue
                        elif i.end > inter.end and i.begin <= inter.end:
                            if i.end-inter.end >= (i.end-i.begin)/2:
                                continue
                            else:
                                all_cont += 1
                new_inter = Interval(inter.begin,inter.end,all_cont)
                wintree[key].remove(inter)
                wintree[key].add(new_inter)
    return wintree


def statGeneID(chr_tree,wintree):
    for key in wintree.keys():
        if 'ctg' not in key and 'tig' not in key and 'utg' not in key:
            for inter in sorted(wintree[key]):
                all_cont = []
                if chr_tree[key].overlap(inter):

                    for i in chr_tree[key].overlap(inter):
                        if i.begin >= inter.begin and i.end < inter.end:
                            all_cont.append(i.data)
                        elif i.end >= inter.begin and i.begin < inter.begin:
                            if i.end-inter.begin >= (i.end-i.begin)/2:
                                all_cont.append(i.data)
                            else:
                                continue
                        elif i.end > inter.end and i.begin <= inter.end:
                            if i.end-inter.end >= (i.end-i.begin)/2:
                                continue
                            else:
                                all_cont.append(i.data)
                new_inter = Interval(inter.begin,inter.end,all_cont)
                wintree[key].remove(inter)
                wintree[key].add(new_inter)
    return wintree


def statGeneLen(chr_tree,wintree):
    for key in wintree.keys():
        if 'ctg' not in key and 'tig' not in key and 'utg' not in key:
            for inter in sorted(wintree[key]):
                all_cont = 0
                if chr_tree[key].overlap(inter):
                    for i in chr_tree[key].overlap(inter):
                        if i.begin >= inter.begin and i.end < inter.end:
                            all_cont += i.data
                        elif i.end >= inter.begin and i.begin < inter.begin:
                            if i.end-inter.begin >= (i.end-i.begin)/2:
                                all_cont += i.data

                            else:
                                continue
                        elif i.end > inter.end and i.begin <= inter.end:
                            if i.end-inter.end >= (i.end-i.begin)/2:
                                continue
                            else:
                                all_cont += i.data
                new_inter = Interval(inter.begin,inter.end,all_cont)
                wintree[key].remove(inter)
                wintree[key].add(new_inter)
    return wintree


def statGCcount(fastaD,wintree):
    for chr in wintree.keys():
        for inter in sorted(wintree[chr]):
            GCnum = 0
            seq = str(fastaD[chr])[inter.begin:inter.end]
            for i in seq:
                if i == 'G' or i == 'C':
                    GCnum += 1
            new_inter = Interval(inter.begin, inter.end, round(GCnum/(inter.end-inter.begin),2))
            wintree[chr].remove(inter)
            wintree[chr].add(new_inter)

    return wintree


def output_result(newtree_D,outputfile=None):
    for key in sorted(newtree_D.keys()):
        if 'ctg' not in key and 'tig' not in key and 'utg' not in key:
            for inter in sorted(newtree_D[key]):
                print(key,inter.begin,inter.end,inter.data,sep='\t',file=outputfile)


def output_result2(newtree_D,outputfile=None):
    for key in sorted(newtree_D.keys()):
        if 'ctg' not in key and 'tig' not in key and 'utg' not in key:
            for inter in sorted(newtree_D[key]):
                print(key,inter.begin,inter.end,','.join(inter.data),sep='\t',file=outputfile)


#outside command 
def MakeWindows(args):
    '''
    Create window files for sliding window analysis
    >>> %(prog)s -f <fasta> -w windows_size -s step_size -o output [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=MakeWindows.__name__,
                        description=MakeWindows.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')


    pReq.add_argument('-f','--fasta',required=True,
            help='Input fasta file')
    pReq.add_argument('-w','--window',required=True,
            type=int,help='Specify window size')
    pReq.add_argument('-s','--step',required=True,
            type=int,help='Specify step length')
    pReq.add_argument('-o','--output',required=True,
            help='Specify the output file name')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.fasta)
    seqD = Fa2dict(args.fasta)
    seq2len = stat_seqLen(seqD)
    sliding_windows(seq2len,args.window,args.step,args.output)
    log.info('Completed! the output file is `{}`'.format(args.output))


def SlidingWindow2GeneNum(args):
    '''
    Counting the number of genes on a sequence using sliding window analysis
    >>> %(prog)s -g <gff_file> -w windows_file [-o output] [Options]
    Note: The gff file must have 'ID=' in the attributes
    '''
    install()
    p = argparse.ArgumentParser(prog=SlidingWindow2GeneNum.__name__,
                        description=SlidingWindow2GeneNum.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('-g','--gff',required=True,
            help='Input gff3 file')
    pReq.add_argument('-t','--type',required=True,
            help='Input the type of extracted gene',
            choices=['ID','Name'])
    pReq.add_argument('-w','--window',required=True,
            help='Input window file')    
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            # default=sys.stdout, 
            help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.gff)
    chr2tree = gff2tree_id(args.gff,args.type)
    check_file_exists(args.window)
    win2tree = window2tree(args.window)
    newtree = statGeneNum(chr2tree,win2tree)
    
    if args.output:
        output_result(newtree,args.output)
        log.info('Completed! the output file is `{}`'.format(args.output.name))
    else:
        output_result(newtree)


def SlidingWindow2GeneID(args):
    '''
    Stating the ID of genes on a sequence using sliding window analysis
    >>> %(prog)s -g <gff_file> -w windows_file [-o output] [Options]
    Note: The gff file must have 'ID=' in the attributes !
    '''
    install()
    p = argparse.ArgumentParser(prog=SlidingWindow2GeneNum.__name__,
                        description=SlidingWindow2GeneNum.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('-g','--gff',required=True,
            help='Input gff3 file')
    pReq.add_argument('-t','--type',required=True,
            help='Input the type of extracted gene',
            choices=['ID','Name'])
    pReq.add_argument('-w','--window',required=True,
            help='Input window file')    
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            # default=sys.stdout, 
            help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.gff)
    chr2tree = gff2tree_id(args.gff,args.type)
    check_file_exists(args.window)
    win2tree = window2tree(args.window)
    newtree = statGeneID(chr2tree,win2tree)
    
    if args.output:
        output_result2(newtree,args.output)
        log.info('Completed! the output file is `{}`'.format(args.output.name))
    else:
        output_result(newtree)


def SlidingWindow2GeneLen(args):
    '''
    Counting the length of genes on a sequence using sliding window analysis
    >>> %(prog)s -g <gff_file> -w windows_file [-o output] [Options]
    Note: The gff file must have 'ID=' in the attributes
    '''
    install()
    p = argparse.ArgumentParser(prog=SlidingWindow2GeneLen.__name__,
                        description=SlidingWindow2GeneLen.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('-g','--gff',required=True,
            help='Input gff3 file')
    pReq.add_argument('-w','--window',required=True,
            help='Input window file')    
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            # default=sys.stdout, 
            help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.gff)
    chr2tree = gff2tree_len(args.gff)
    check_file_exists(args.window)
    win2tree = window2tree(args.window)
    newtree = statGeneLen(chr2tree,win2tree)
    if args.output:
        output_result(newtree,args.output)
        log.info('Completed! the output file is `{}`'.format(args.output.name))
    else:
        output_result(newtree)


def SlidingWindow2GCcontent(args):
    '''
    Counting GC content on a sequence using sliding window analysis
    >>> %(prog)s -f <fasta> -w windows_file [-o output] [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=SlidingWindow2GCcontent.__name__,
                        description=SlidingWindow2GCcontent.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('-f','--fasta',required=True,
            help='Input fasta file')
    pReq.add_argument('-w','--window',required=True,
            help='Input window file')    
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            # default=sys.stdout, 
            help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.fasta)
    fastaD = Fa2dict(args.fasta)
    check_file_exists(args.window)
    win2tree = window2tree(args.window)
    newtree = statGCcount(fastaD,win2tree)
    
    if args.output:
        output_result(newtree,args.output)
        log.info('Completed! the output file is `{}`'.format(args.output.name))
    else:
        output_result(newtree)


def main():
    actions = (
            ("MakeWindows", "Create window files for sliding window analysis"),
            ("SlidingWindow2GeneNum", "Counting the number of genes on a sequence using sliding window analysis"),
            ("SlidingWindow2GeneID", "Stating the ID of genes on a sequence using sliding window analysis"),
            ("SlidingWindow2GeneLen", "Stating the length of genes on a sequence using sliding window analysis"),
            ("SlidingWindow2GCcontent", "Stating GC content on a sequence using sliding window analysis"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
