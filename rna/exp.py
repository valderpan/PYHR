#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/05/18


import re
import argparse
import pandas as pd
from path import Path
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, check_path_exists, richlog
from PYHR.apps.base import listify


log = richlog()

def getFpkmfile(file_path):
    files = [file for file in Path(file_path).files() if file.endswith('.tab')]
    return files


def merge2matrix(files,valuetype,genetype,output_file):
    mergedf_list = []
    for file in files:
        file_name = re.findall('([0-9A-Za-z\_\.\-]*)\.tab', file)[0]
        check_file_exists(file)
        df = pd.read_table(file,sep='\t',float_precision='round_trip')
        if genetype == 'ID':
            df2merge = df.loc[:,['Gene ID',valuetype]]
            df2merge['value'] = round(df2merge[valuetype],2)
            df2merge = df2merge.loc[:,['Gene ID','value']]
            df2merge.rename(columns={'value':file_name},inplace=True)
            mergedf_list.append(df2merge)
        elif genetype == 'Name':
            df2merge = df.loc[:,['Gene Name',valuetype]]
            df2merge['value'] = round(df2merge[valuetype],2)
            df2merge = df2merge.loc[:,['Gene Name','value']]
            df2merge.rename(columns={'value':file_name},inplace=True)
            mergedf_list.append(df2merge)
        else:
            log.error(f'Option {genetype} does not exist, it can only be one of [ID or Name]')
    if genetype == 'ID':
        res = mergedf_list[0].sort_values(by='Gene ID')
        for i in range(1,len(mergedf_list)):
            res = pd.merge(res,mergedf_list[i],on="Gene ID")
    elif genetype == 'Name':
        res = mergedf_list[0].sort_values(by='Gene Name')
        for i in range(1,len(mergedf_list)):
            res = pd.merge(res,mergedf_list[i],on="Gene Name")
    else:
        log.error(f'Option {genetype} does not exist, it can only be one of [ID or Name]')
    # for i in range(1,len(mergedf_list)):
    #     res = pd.merge(res,mergedf_list[i],on="Gene ID")
    log.info('The results file is output to `{}` '.format(output_file))
    res.to_excel(output_file,header=True,index=False)


#outside command
def stringtie2ExpMatrix(args):
    '''
    Combine the expression values of each sample into a matrix after calculating the expressions from Hisat2+stringtie pipeline
    >>> %(prog)s <Path of tab files> <FPKM|TPM> -o <output.xlsx>
    '''
    install()
    p = argparse.ArgumentParser(prog=stringtie2ExpMatrix.__name__,
                        description=stringtie2ExpMatrix.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('path', 
            help='Path to store stringtie outpur file (The target file under the path must end with `.tab`)')
    pReq.add_argument('-v','--valuetype',choices=['FPKM','TPM'] ,
            help='Select the type of data to be extracted')
    pReq.add_argument('-g','--genetype',choices=['ID','Name'] ,
            help='Select the type of data to be extracted,{ID}:ENSEMBL,{Name}:SYMBL')
    pReq.add_argument('-o', '--output', required=True,
            help='Specify the output file (.xlsx)')
    pOpt.add_argument('-h', '--help', action='help',
            help='Show help message and exit.')

    args = p.parse_args(args)

    check_file_exists(args.path)
    files = getFpkmfile(args.path)
    log.info('Start merge ...')
    merge2matrix(files,args.valuetype,args.genetype,args.output)


def main():
    actions = (
            ("stringtie2ExpMatrix", "Combine the expression values of each sample into a matrix after calculating the expressions from Hisat2+stringtie pipeline"),
            #("concatKs2ggplotdensity", "Extract query seq by header"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()