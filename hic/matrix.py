#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2023/01/01


import re
import sys
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

def cleanmatrix(bed,matrix,chrom):
    beddf = pd.read_table(bed,names=['chrom','start','end','bnum'],converters={'chrom':str})
    matrixdf = read_file(matrix,['b1','b2','contact'])
    if chrom not in [str(x) for x in beddf['chrom'].unique().tolist()]:
        log.error(f'{chrom} does not exist,Please check the input chrom parameter!')
        sys.exit()
    else:
        qbeddf = beddf[beddf['chrom'] == chrom].drop_duplicates().reset_index(drop=True)
        minN = qbeddf['bnum'].min()
        maxN = qbeddf['bnum'].max()
        N = maxN - minN + 1
        qmatrixdf = matrixdf[(matrixdf['b1'] >= minN) & (matrixdf['b1'] < maxN) & (matrixdf['b2'] >= minN) & (matrixdf['b2'] < maxN)]
        qmatrixdf['b1'] = qmatrixdf['b1'] - minN
        qmatrixdf['b2'] = qmatrixdf['b2'] - minN
        qmatrixdf.index = np.array(np.arange(qmatrixdf.shape[0]))
        sparse_matrix= sparse.coo_matrix((qmatrixdf['contact'], (qmatrixdf['b1'], qmatrixdf['b2'])), shape=(N,N), dtype=float).toarray()
    return sparse_matrix,N,qbeddf


def asymmetric_matrix(sparse_matrix):
    diag_matrix   = np.diag(np.diag(sparse_matrix))
    as_sparse_matrix = sparse_matrix.T + sparse_matrix 
    as_sparse_matrix = as_sparse_matrix-diag_matrix
    return as_sparse_matrix


def createNNmatrix(sparse_matrix,genome,N,beddf,outdir,sample,chrom,resolution):
    sparse_matrix = pd.DataFrame(sparse_matrix)
    sparse_matrix = sparse_matrix.fillna(0)
    inddf         = np.arange(N)
    headers_ref   = [genome for x in inddf]
    bin_num_df    = pd.Series(beddf['bnum']).apply(lambda x:str(x))
    headers_ref   = pd.Series(headers_ref)
    chromdf       = pd.Series(beddf['chrom'].apply(lambda x:str(x)))
    startdf       = pd.Series(beddf['start']).apply(lambda x:str(x))
    enddf         = pd.Series(beddf['end']).apply(lambda x:str(x))
    headers_suf   = chromdf.str.cat(startdf,sep=':')
    headers_suf   = headers_suf.str.cat(enddf,sep="-")
    headers       = bin_num_df.str.cat([headers_ref,headers_suf],sep="|")
    headers       = list(headers)
    sparse_matrix.columns=headers
    sparse_matrix.index=headers
    output = f"{outdir}/{sample}.{genome}_{chrom}.{resolution}.matrix"
    sparse_matrix.to_csv(f'{output}',sep="\t",index=True,header=True,float_format='%0.3f')
    return output

def buildmatrix_basesize(dirpath,bsize,fasize):
    hicdir = [d for d in Path(dirpath).dirs() if d.basename() == 'hic_results'][0]
    for d in Path(hicdir).dirs():
        if d.basename() == 'data':
            sampleID = [i for i in Path(d).dirs()]
            tmp = ','.join([i.basename() for i in sampleID])
            log.info(f'A total of {len(sampleID)} samples were identified as [{tmp}]')
            for s in sampleID:
                validfile = [file for file in Path(s).files() if file.endswith('.allValidPairs')][0]
                cmd1 = f'mkdir -p {hicdir}/matrix/{s.basename()}/iced/{bsize}'
                cmd2 = f'mkdir -p {hicdir}/matrix/{s.basename()}/raw/{bsize}'
                cmd3_1 = f"cat {validfile} | build_matrix --matrix-format upper "
                cmd3_2 = f"--binsize {bsize} --chrsizes {fasize} --ifile /dev/stdin "
                cmd3_3 = f"--oprefix {hicdir}/matrix/{s.basename()}/raw/{bsize}/{s.basename()}_{bsize}"
                cmd3 = f"{cmd3_1}{cmd3_2}{cmd3_3}"
                
                cmd4_1 = f"ice --results_filename {hicdir}/matrix/{s.basename()}/iced/{bsize}/{s.basename()}_{bsize}_iced.matrix "
                cmd4_2 = "--filter_low_counts_perc 0.02 --filter_high_counts_perc 0 --max_iter 100 --eps 0.1 --remove-all-zeros-loci --output-bias 1 "
                cmd4_3 = f"--verbose 1 {hicdir}/matrix/{s.basename()}/raw/{bsize}/{s.basename()}_{bsize}.matrix"
                cmd4 = f"{cmd4_1}{cmd4_2}{cmd4_3}"
                runshell(cmd1)
                runshell(cmd2)
                runshell(cmd3)
                runshell(cmd4)


def icedmatrix(ice,matrix):
    prefix = re.findall('([0-9a-zA-Z\_\-\.\/]+)\.matrix',matrix)
    if prefix != 0:
        cmd = f"python {ice} --results_filename {prefix[0]}_iced.matrix \
            --filter_low_counts_perc 0.02 --filter_high_counts_perc 0 \
            --max_iter 100 --eps 0.1 --remove-all-zeros-loci --output-bias 1 \
            {matrix}"
        runshell(cmd)
    else:
        log.error(f'File prefix extraction error, the obtained prefix was {prefix}')
    return f"{prefix[0]}_iced.matrix"

## outside command
def ICEMatrix(args):
    """
    Use ice to process the raw matrix of HiC-Pro
    >>> %(prog)s <matrix> [iced.matrix][Options]
    """ 
    install()
    p = argparse.ArgumentParser(prog=ICEMatrix.__name__,
                        description=ICEMatrix.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('ice',
            type=str,
            help='Input the path of ice')
    pReq.add_argument('matrix',
            type=str,
            help='Input HiC-Pro cisSparse raw matrix file')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.matrix)
    outfile = icedmatrix(args.ice,args.matrix)
    log.info(f'The results file is output to `{outfile}`')


def CreateNNMatrix(args):
    """
    Create dense matrix from HiC-Pro results
    Referenced from Condense_CisMatrix.py and https://github.com/nservant/HiC-Pro/blob/master/bin/utils/sparseToDense.py
    >>> %(prog)s -b <bed> -m <matrix> -g genome_name -c chrom_name -o outdir -s sample_name -r resolution [Options]
    """ 
    install()
    p = argparse.ArgumentParser(prog=CreateNNMatrix.__name__,
                        description=CreateNNMatrix.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-b','--bed',
            dest='bed',type=str,required=True,
            help='Input HiC-Pro bed file')
    pReq.add_argument('-m','--matrix',
            dest='matrix',type=str,required=True,
            help='Input HiC-Pro cisSparse matrix file')
    pReq.add_argument('-g','--genome',
            dest='genome',type=str,required=True,
            help='Specify the name of the genome to be used,such as hg19,hg18...')
    pReq.add_argument('-c','--chrom',
            dest='chrom',type=str,required=True,
            help='Specify the chromosome ID')
    pReq.add_argument('-o','--outdir',
            dest='outdir',type=str,required=True,
            help='Specify the output path')
    pReq.add_argument('-s','--sample',
            dest='sample',type=str,required=True,
            help='Specify the sample ID')
    pReq.add_argument('-r','--resolution',
            dest='resolution',type=str,required=True,
            help='Specify the resolution size')
    pOpt.add_argument('-asym', '--asymmetric',
            dest='asymmetric',default="no",type=str,choices=["yes","no"],
            help='if need symmetric yes else no')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)

    check_file_exists(args.bed)
    check_file_exists(args.matrix)
    sparse_matrix,N,qbeddf = cleanmatrix(args.bed,args.matrix,args.chrom)
    if args.asymmetric == "yes":
        sparse_matrix = asymmetric_matrix(sparse_matrix)
    outfile = createNNmatrix(sparse_matrix,args.genome,N,qbeddf,args.outdir,args.sample,args.chrom,args.resolution)
    log.info(f'The results file is output to `{outfile}`')


def bulidMatrix(args):
    """
    Create the HiC-Pro result matrix based on the given resolution
    >>> %(prog)s -d <HiC-Pro result dir> -s binsize -f fasta.size [Options]
    """ 
    install()
    p = argparse.ArgumentParser(prog=bulidMatrix.__name__,
                        description=bulidMatrix.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-d','--dirpath',
            dest='dirpath',type=str,required=True,
            help='Input HiC-Pro result dir')
    pReq.add_argument('-s','--binsize',
            dest='binsize',type=int,required=True,
            help='Input the bin size, also called resolution')
    pReq.add_argument('-f','--fastasize',
            dest='fastasize',type=str,required=True,
            help='Input the fasta.size used for HiC-Pro run.')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    args = p.parse_args(args)
    
    check_file_exists(args.dirpath)
    check_file_exists(args.fastasize)
    buildmatrix_basesize(args.dirpath,args.binsize,args.fastasize)


def main():
    actions = (
            ("ICEMatrix", "Use ice to process the raw matrix of HiC-Pro"),
            ("CreateNNMatrix", "Create dense matrix from HiC-Pro results"),
            ("bulidMatrix", "Create the HiC-Pro result matrix based on the given resolution")
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
