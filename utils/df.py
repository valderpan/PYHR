#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/05/18


import sys
import argparse
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify, read_file


log = richlog()

def read2List(file):
    '''
    Store all values from the file in the list(values are split by "\n")
    '''
    with open(file) as f:
        _vL = [l for l in (line.strip() for line in f)]
    return _vL

def subdf(vL,df,col_name,col_type):
    if col_type == 'name':
        if col_name not in [col for col in df]:
            log.error('{} is not in the column name of dataframe ! Please check if the input ab is correct !'.format(col_name))
            sys.exit()
        else:
            sdf = df[df[col_name].isin(vL)]
    elif col_type == 'index':
        if int(col_name) > len(df.columns.values.tolist()):
            log.error('Index {} is not in the column name of dataframe ! Please check if the input ab is correct !'.format(col_name))
            sys.exit()
        else:
            sdf = df[df.iloc[:,int(col_name)].isin(vL)]
    return sdf


## outside command 
def extractSubdf(args):
    """
    Extract the sub dataframe based on the values inside the list
    >>> %(prog)s <valuefile> <dataframe> <colname> -o <output> [Options]
    """ 
    install()
    p = argparse.ArgumentParser(prog=extractSubdf.__name__,
                        description=extractSubdf.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('valuefile', 
            help='Input file with values stored by `\\n` partition')
    pReq.add_argument('dataframe', 
            help='Input the dataframe file to be extracted\nNote: This file must have a header !!!')
    pReq.add_argument('-c','--colname', required=True,
            help='Specify which column of the dataframe file to extract based on')
    pReq.add_argument('-t','--coltype', required=True,choices=['name','index'],
            help='Specify column type of the dataframe file to extract based on, \
                it can be a column name or a column index\nNote:Use only if the -c parameter is present')
    pReq.add_argument('-o', '--output', required=True,
            help='The name of the output file')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    check_file_exists(args.valuefile)
    vl = read2List(args.valuefile)
    check_file_exists(args.dataframe)
    df = read_file(args.dataframe)
    res = subdf(vl,df,args.colname,args.coltype)
    log.info('The results file is output to {}'.format(args.output))
    res.to_csv(args.output,header=False,index=False,sep='\t')


def main():
    actions = (
            ("extractSubdf", "Extract the sub dataframe based on the values inside the list"),
            #("concatKs2ggplotdensity", "Extract query seq by header"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()