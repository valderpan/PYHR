#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/05/19


import re
import sys
import argparse
import pandas as pd
from path import Path
from natsort import natsorted
from collections import OrderedDict
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify, read_file
from PYHR.utils.fasta import Fa2dict


log = richlog()

def read_bed(bed):
    pre = re.findall('([0-9a-zA-Z\-\_\.]+)\.bed',bed)[0]
    df = pd.read_table(bed,sep='\t',names=['chr','start','end','index'],float_precision='round_trip')
    return df,pre


def read_matrix(matrix):
    pre = re.findall('([0-9a-zA-Z\-\_\.]+)\.matrix',matrix)[0]
    df = pd.read_table(matrix,sep='\t',names=['bin1','bin2','score'],float_precision='round_trip')
    return df,pre


def modify_icedmatrix(matrixdf,pre,output=None):
    matrixdf['newbin1'] = matrixdf['bin1'].astype(int)+1
    matrixdf['newbin2'] = matrixdf['bin2'].astype(int)+1
    newdf = matrixdf.loc[:,['newbin1','newbin2','score']]
    if output:
        newdf.to_csv(output,sep='\t',header=False,index=False)
        log.info('The results file is output to `{}` '.format(output))
    else:
        newdf.to_csv('{}.modified.matrix'.format(pre),sep='\t',header=False,index=False)
        log.info('The results file is output to `{}.modified.matrix` '.format(pre))


def remove_ctg_matrix(matrixdf,pre,rownum,output=None):
    log.info('The input index of the last bin of the chromosome is {}'.format(rownum))
    matrixdf = matrixdf[matrixdf['bin1'] <= int(rownum)]
    matrixdf = matrixdf[matrixdf['bin2'] <= int(rownum)]
    

    if output:
        matrixdf.to_csv(output,sep='\t',header=False,index=False)
        log.info('The results file is output to `{}` '.format(output))
    else:
        matrixdf.to_csv('{}.noctg.matrix'.format(pre),sep='\t',header=False,index=False)
        log.info('The results file is output to `{}.noctg.matrix` '.format(pre))


# TODO 关于提取SubMatrix这一点还没想好目的是什么
# def get_Qregion(bed_df,chrID,pre):
#     q_df = bed_df[bed_df['chr'] == chrID]
#     i_l = q_df['index'].tolist()
#     i_s = min(i_l)
#     i_e = max(i_l)
#     q_df.to_csv('{}.{}.bed'.format(pre,chrID),sep='\t',header=False,index=False)
#     return i_s,i_e


# def get_Qmatrix(matrix_df,i_s,i_e):
#     q_df = matrix_df[(matrix_df['bin1'].map(int) >= int(i_s)) & (matrix_df['bin1'].map(int) <= int(i_e))]
#     q_df = q_df[(q_df['bin2'].map(int) >= int(i_s)) & (q_df['bin2'].map(int) <= int(i_e))]
#     return q_df


# def output(q_df,chrID,pre):
#     q_df['score'] = round(q_df['score'],6)
#     q_df.to_csv('{}.{}.matrix'.format(pre,chrID),sep='\t',header=False,index=False)


def findresultfile(reslut_path):
    resultdir = [dir for dir in Path(reslut_path).dirs() if 'hic_results' in dir][0]
    stats_dir  = [dir for dir in Path(resultdir).dirs() if 'stats' in dir][0]
    sample_dirs = [dir for dir in Path(stats_dir).dirs()]
    return sample_dirs


def statsresult(dirlist,output):
    from rich import print
    for sample in dirlist:
        item2num = OrderedDict()
        files = Path(sample).files()
        for file in files:
            if file.endswith('mpairstat'):
                with open(file,'r') as f1:
                    lines = (line.strip() for line in f1)
                    for line in lines:
                        if line.startswith('#'):
                            continue
                        else:
                            line_list = line.split('\t')
                            if line_list[0] == 'Total_pairs_processed':
                                item2num['Clean Paired-end Reads'] = line_list[1]
                            elif line_list[0] == 'Unmapped_pairs':
                                item2num['Unmapped Paired-end Reads'] = line_list[1]
                                #item2num['Unmapped Paired-end Reads Rate (%)'] = line_list[2]
                                item2num['Unmapped Paired-end Reads Rate (%)'] = round(int(item2num['Unmapped Paired-end Reads'])/int(item2num['Clean Paired-end Reads'])*100,3)
                            elif line_list[0] == 'Pairs_with_singleton':
                                item2num['Paired-end Reads with Singleton'] = line_list[1]
                                #item2num['Paired-end Reads with Singleton Rate(%)'] = line_list[2]
                                item2num['Paired-end Reads with Singleton Rate(%)'] = round(int(item2num['Paired-end Reads with Singleton'])/int(item2num['Clean Paired-end Reads'])*100,3)
                            elif line_list[0] == 'Multiple_pairs_alignments':
                                item2num['Multi Mapped Paired-end Reads'] = line_list[1]
                                #item2num['Multi Mapped Ratio (%)'] = line_list[2]
                                item2num['Multi Mapped Ratio (%)'] = round(int(item2num['Multi Mapped Paired-end Reads'])/int(item2num['Clean Paired-end Reads'])*100,3)
                            elif line_list[0] == 'Unique_paired_alignments':
                                item2num['Unique Mapped Paired-end Reads'] = line_list[1]
                                #item2num['Unique Mapped Ratio (%)'] = line_list[2]
                                item2num['Unique Mapped Ratio (%)'] = round(int(item2num['Unique Mapped Paired-end Reads'])/int(item2num['Clean Paired-end Reads'])*100,3)
            elif file.endswith('.mRSstat'):
                with open(file,'r') as f2:
                    lines = (line.strip() for line in f2)
                    for line in lines:
                        if line.startswith('#'):
                            continue
                        else:
                            line_list = line.split('\t')
                            if line_list[0] == 'Dangling_end_pairs':
                                item2num['Dangling End Paired-end Reads'] = line_list[1]
                            elif line_list[0] == 'Self_Cycle_pairs':
                                item2num['Self Circle Paired-end Reads'] = line_list[1]
                            elif line_list[0] == 'Dumped_pairs':
                                item2num['Dumped Paired-end Reads'] = line_list[1]
                            elif line_list[0] == 'Valid_interaction_pairs':
                                item2num['Interaction Paired-end Reads'] = line_list[1]
            elif file.endswith('.mergestat'):
                with open(file,'r') as f3:
                    lines = (line.strip() for line in f3)
                    for line in lines:
                        if line.startswith('#'):
                            continue
                        else:
                            line_list = line.split('\t')
                            if line_list[0] == 'valid_interaction_rmdup':
                                item2num['Lib Valid Paired-end Reads'] = line_list[1]
        
        print(sample+':',file=output)
        print('Statistics of mapping',file=output)
        print('\t'+'Clean Paired-end Reads',item2num['Clean Paired-end Reads'],file=output)
        print('\t'+'Unmapped Paired-end Reads',item2num['Unmapped Paired-end Reads'],file=output)
        print('\t'+'Unmapped Paired-end Reads Rate (%)',item2num['Unmapped Paired-end Reads Rate (%)'],file=output)
        print('\t'+'Paired-end Reads with Singleton',item2num['Paired-end Reads with Singleton'],file=output)
        print('\t'+'Paired-end Reads with Singleton Rate(%)',item2num['Paired-end Reads with Singleton Rate(%)'],file=output)
        print('\t'+'Multi Mapped Paired-end Reads',item2num['Multi Mapped Paired-end Reads'],file=output)
        print('\t'+'Multi Mapped Ratio (%)',item2num['Multi Mapped Ratio (%)'],file=output)
        print('\t'+'Unique Mapped Paired-end Reads',item2num['Unique Mapped Paired-end Reads'],file=output)
        print('\t'+'Unique Mapped Ratio (%)',item2num['Unique Mapped Ratio (%)'],file=output)

        print('Statistics of valid reads',file=output)
        print('\t'+'Unique Mapped Paired-end Reads',item2num['Unique Mapped Paired-end Reads'],file=output)
        print('\t'+'Dangling End Paired-end Reads',item2num['Dangling End Paired-end Reads'],file=output)
        print('\t'+ 'Dangling End Rate (%)',round(int(item2num['Dangling End Paired-end Reads']) / int(item2num['Unique Mapped Paired-end Reads']) *100, 3),file=output)
        print('\t'+ 'Self Circle Paired-end Reads', item2num['Self Circle Paired-end Reads'],file=output)
        print('\t'+ 'Self Circle Rate (%)',round(int(item2num['Self Circle Paired-end Reads']) / int(item2num['Unique Mapped Paired-end Reads']) *100, 3),file=output)
        print('\t'+ 'Dumped Paired-end Reads', item2num['Dumped Paired-end Reads'],file=output)
        print('\t'+ 'Dumped Rate (%)', round(int(item2num['Dumped Paired-end Reads']) / int(item2num['Unique Mapped Paired-end Reads']) *100, 3),file=output)
        print('\t'+ 'Interaction Paired-end Reads', item2num['Interaction Paired-end Reads'],file=output)
        print('\t'+ 'Interaction Rate (%)',round(int(item2num['Interaction Paired-end Reads']) / int(item2num['Unique Mapped Paired-end Reads']) *100 , 3),file=output)
        print('\t'+ 'Lib Valid Paired-end Reads', item2num['Lib Valid Paired-end Reads'],file=output)
        Lib_Valid_Rate = round(int(item2num['Lib Valid Paired-end Reads']) / int(item2num['Interaction Paired-end Reads']) *100, 3)
        print('\t', 'Lib Valid Rate (%)',Lib_Valid_Rate,file=output)
        print('\t','Lib Dup (%)',round(100-Lib_Valid_Rate,3),file=output)

#TODO
def addmissingconcat(beddf,matrix,output):
    bD1 = {};approw = []
    for _,row in matrix.iterrows():
        bD1[str(int(row['bin1']))+'_'+str(int(row['bin2']))] = row['score']
    for key in bD1.keys():
        rowdf = pd.DataFrame({"bin1":key.split('_')[1],"bin2":key.split('_')[0],"score":bD1[key]})
        approw.append(rowdf)
    ndf = pd.concat(approw)
    resdf = pd.concat(matrix,ndf,axis=1)
    resdf = resdf.sort_values(by=['bin1','bin2'])
    resdf.to_csv(output,sep='\t',index=False,header=True)

## outside command 
def modifyMatrixIndex(args):
    """
    The bin index in HiC-Pro iced_matirx should start from 1 instead of 0
    >>> %(prog)s <matirx> [Options] -o output
    """
    install()
    p = argparse.ArgumentParser(prog=modifyMatrixIndex.__name__,
                        description=modifyMatrixIndex.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('matrix', 
            help='Input the matrix file ')
    pOpt.add_argument('-o', '--output',default=None,
            help='output file [default : prefix.modified.matrix]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    check_file_exists(args.matrix)
    matrixdf,pre = read_matrix(args.matrix)
    modify_icedmatrix(matrixdf,pre,args.output)


def removeCtgMatrix(args):
    """
    Delete the contig level interaction signal from the Hi-C matrix file
    ==============================================================================
    row num : the index value corresponding to the last bin of the last chromosome
    ==============================================================================
    >>> %(prog)s <matirx> <rownum> [Options] -o output
    """
    install()
    p = argparse.ArgumentParser(prog=removeCtgMatrix.__name__,
                        description=removeCtgMatrix.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('matrix', 
            help='Input the matrix file ')
    pReq.add_argument('rownum', type=int,
            help='Input the index of the last bin of the chromosome')
    pOpt.add_argument('-o', '--output',default=None,
            help='output file  [default : prefix.noctg.matrix]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    check_file_exists(args.matrix)
    matrixdf,pre = read_matrix(args.matrix)
    remove_ctg_matrix(matrixdf,pre,args.rownum,args.output)


# TODO
# def extractSubMatrix(args):
#     """
#     Delete the contig level interaction signal from the Hi-C matrix file
#     row num : the index value corresponding to the last bin of the last chromosome
#     >>> %(prog)s <matirx> <rownum> [Options] -o output
#     """
#     install()
#     p = argparse.ArgumentParser(prog=extractSubMatrix.__name__,
#                         description=extractSubMatrix.__doc__,
#                         formatter_class=argparse.RawTextHelpFormatter,
#                         conflict_handler='resolve')
#     pReq = p.add_argument_group('Required arguments')
#     pOpt = p.add_argument_group('Optional arguments')
#     pReq.add_argument('matrix', 
#             help='Input the matrix file ')
#     pReq.add_argument('rownum', type=int,
#             help='Input the index of the last bin of the chromosome')
#     pOpt.add_argument('-o', '--output',default=None,
#             help='output file  [default : prefix.noctg.matrix]')
#     pOpt.add_argument('-h', '--help', action='help',
#             help='show help message and exit.')
    
#     args = p.parse_args(args)


def statHiCPro(args):
    """
    Convert HiC-Pro output stats results into publishable level tables\n
    ====================================================================
    Note: This script is available in both absolute and relative paths !
    ====================================================================
    >>> %(prog)s <path> [Options] -o output
    """
    install()
    p = argparse.ArgumentParser(prog=statHiCPro.__name__,
                        description=statHiCPro.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    
    pReq.add_argument('path', 
            help='Input the HiC-Pro results path [Directory obtained by HiC-Pro -o parameter] ')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    check_file_exists(args.path)
    sample_dirs = findresultfile(args.path)
    if len(sample_dirs) == 0:
        log.error('The contents of the {} directory are incorrect and the correct file cannot be found'.format(args.path))
        sys.exit()
    else:
        statsresult(sample_dirs,args.output)

def main():
    actions = (
            ("modifyMatrixIndex", "The bin index in HiC-Pro iced_matirx should start from 1 instead of 0"),
            ("removeCtgMatrix", "Delete the contig level interaction signal from the Hi-C matrix file"),
            ("statHiCPro", "Convert HiC-Pro output stats results into publishable level tables\n"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()