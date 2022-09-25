#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/05/17


import re
import sys
import math
import argparse
import pandas as pd
from path import Path
from natsort import natsorted
from rich.traceback import install
from intervaltree import IntervalTree,Interval
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify


log = richlog()

def cal_Ks_valuecounts(Ks_file,limit):
    df = pd.read_table(Ks_file,sep='\t')
    ksdf = df.iloc[:,[0,3]]
    if 'dS-ng' not in [column for column in df]:
        log.error('Please enter the correct KaKs result file !(The correct file should be the result obtained from the post-processing of the synonymous_calc.py calculation)')
        sys.exit()
    else:
        ksdf = ksdf[(ksdf['dS-ng'] >= 0) & (ksdf['dS-ng'] < limit)]
        ksdf.iloc[:,1] = round(ksdf.iloc[:,1],3)
        ks_count_dict = ksdf.iloc[:, 1].value_counts().to_dict()
        return ks_count_dict


def split_windows(L,S):
    limit = float(L)
    step = float(S)
    tree = IntervalTree()
    windows = math.ceil(limit / step)
    win_s = 0
    for i in range(windows):
        if i == 0:
            win_e = win_s + step
            tree.addi(win_s, win_e, 0)
        else:
            if i < windows - 1:
                win_s += step
                if win_s + step < limit:
                    win_e = win_s + step
                else:
                    win_e = limit
                tree.addi(win_s, win_e, 0)
            else:
                win_s += step
                win_e = limit
                tree.addi(win_s, win_e, 0)
    return tree


def cal_Ks_density(ks_count_dict,tree,limit):
    for key in ks_count_dict.keys():
        inter_L = sorted(tree[key])
        if len(inter_L) == 1:
            inter = inter_L[0]
            new_inter = Interval(inter.begin, inter.end, inter.data + ks_count_dict[key])
            tree.remove(inter)
            tree.add(new_inter)
        else:
            if key == limit:
                inter = sorted(tree)[-1]
                new_inter = Interval(inter.begin,inter.end,inter.data +ks_count_dict[key])
                tree.remove(inter)
                tree.add(new_inter)
            else:
                print('This value lies in two intervals!!!')
    return tree


def convert_result(newtree, ks_name):
    Dict= {}
    for inter in sorted(newtree):
        Dict[round(inter.begin, 2)] = inter.data
    df = pd.DataFrame([Dict]).T
    df.rename(columns={0: '{}_counts'.format(ks_name)}, inplace=True)
    df[ks_name] = df['{}_counts'.format(ks_name)] / sum(df['{}_counts'.format(ks_name)]) * 100
    output = df.sort_index().reset_index()
    return output


def run_Ksdensity(Ks_file,limit,step,ks_name):
    ks_count_dict = cal_Ks_valuecounts(Ks_file,limit)
    tree = split_windows(limit,step)
    newtree = cal_Ks_density(ks_count_dict,tree,limit)
    output = convert_result(newtree,ks_name)
    return output


def read_Ksfile(ks_file_path):
    ks_file_list = [i for i in Path(ks_file_path).files() if i.endswith('KaKs.result')]
    return natsorted(ks_file_list)


def merge_ks_reads_mapping_df(ks_files,Upper_limit,step):
    file_dict = {}
    ks_file_name_list = []
    for i in ks_files:
        j = re.findall('([A-Za-z0-9\_]+)\.KaKs.result',i)[0]
        ks_file_name_list.append(j)
    log.info('Read first Ks file:{}...'.format(ks_files[0]))
    first_file = run_Ksdensity(ks_files[0],Upper_limit,step,ks_file_name_list[0])
    for i in range(1,len(ks_files)):
        log.info('Read Ks file:{}...'.format(ks_files[i]))
        file_dict[i] = run_Ksdensity(ks_files[i],Upper_limit,step,ks_file_name_list[i])
        if i < len(ks_files):
            first_file = pd.merge(first_file,file_dict[i],on='index',how='outer')
    choose_columns_number =  []
    for i in range(0,first_file.shape[1],2):
        choose_columns_number.append(i)
    output_ks_file = first_file.iloc[:,choose_columns_number]
    output_ks_file = output_ks_file.fillna(0)
    output_ks_file = output_ks_file.sort_values(by='index')
    return output_ks_file


def output_results(output_ks_file,output_name):
    log.info('The results file is output to {}'.format(output_name))
    output_ks_file.to_excel(output_name,header=True,index=False)


def Read_files(files):
    C_L = []
    for file in files:
        check_file_exists(file)
        name = re.findall('([0-9a-zA-Z\-\_]+)\.KaKs',file)[0]
        df = pd.read_table(file,sep='\t')
        df['name'] = name
        df = df[df['dS-ng'] > 0]
        qdf = df.iloc[:,[0,3]]
        C_L.append(qdf)
    return  C_L


def ConcatOutput(c_l,output):
    result = pd.concat(c_l)
    log.info('The results file is output to `{}` '.format(output))
    result.to_csv(output,sep='\t',header=['Name','Ks'],index=False)



#outside command 
def calKsdist(args):
    '''
    Calculate the interval distribution of gene pairs with different synonymous substitution ratio
    >>> %(prog)s -p <KaKs_path> -l limit -s step -o <output.xlsx> [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=calKsdist.__name__,
                        description=calKsdist.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')


    pReq.add_argument('-p', '--path', required=True, 
            help='Path to store ks file (The target file under the path must end with `.KaKs.result`)')
    pReq.add_argument('-l', '--limit', required=True, 
            type=float,help='Upper limit of ks reads mapping')
    pReq.add_argument('-s', '--step',required=True,
            type=float,help="Input the step size")
    pReq.add_argument('-o', '--output', required=True,
            help='The name of the output file(.xlsx)')
    pOpt.add_argument('-h', '--help', action='help',
            help='Show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.path)

    ks_files = read_Ksfile(args.path)
    print(ks_files)
    result = merge_ks_reads_mapping_df(ks_files,args.limit,args.step)
    output_results(result,args.output)


def concatKs2ggplotdensity(args):
    '''
    Extract the Ks values from each file, merge them and convert them to ggplot2 friendly format
    >>> %(prog)s <KaKs_path> -o <output.csv/tab>
    '''
    install()
    p = argparse.ArgumentParser(prog=concatKs2ggplotdensity.__name__,
                        description=concatKs2ggplotdensity.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('path', 
            help='Path to store ks file (The target file under the path must end with `.KaKs.result`)')
    pReq.add_argument('-o', '--output', required=True,
            help='The name of the output file(.csv or .tab)')
    pOpt.add_argument('-h', '--help', action='help',
            help='Show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.path)
    ks_files = read_Ksfile(args.path)
    c_l = Read_files(ks_files)
    ConcatOutput(c_l,args.output)


def main():
    actions = (
            ("calKsdist", "Calculate the interval distribution of gene pairs with different synonymous substitution ratio"),
            ("concatKs2ggplotdensity", "Extract query seq by header"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()