#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2025/02/21

import os
import re
import sys
import gzip
import argparse
import pandas as pd
from path import Path
from rich.traceback import install
from rich.console import Console
from rich.table import Column, Table
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify, read_file, runshell

log = richlog()
console = Console()


def read_rename_file(rename_file):
    """
    读取替换文件，返回一个dict对象
    替换文件格式为(以\t分割)：
    old_id1 new_id1
    old_id2 new_id2
    ...
    注意：默认没有header！
    """
    replace_dict = {}
    with open(rename_file, 'r') as f:
        lines = (line.strip() for line in f)
        for line in lines:
            old_id, new_id = line.split('\t')
            replace_dict[old_id] = new_id
    return replace_dict


def rename_vcf_header(vcf_file,replace_dict,output_vcf):
    """
    读取并替换VCF文件中的sample ID
    此函数只处理VCF的#CHROM行，不对其他行进行处理，因此没有额外的读取vcf的function
    """
    with gzip.open(vcf_file, 'rt') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#CHROM'):
                columns = line.strip().split('\t')
                # 找到样本ID的起始列
                sample_start_index = columns.index('FORMAT') + 1
                # 替换样本ID
                for i in range(sample_start_index, len(columns)):
                    old_id = columns[i]
                    if old_id in replace_dict:
                        columns[i] = replace_dict[old_id]
                # 写入替换后的行
                new_line = '\t'.join(columns) + '\n'
                outfile.write(new_line)
            else:
                # 非样本ID行直接写入
                outfile.write(line)


def bgzip_vcf(vcf_file):
    """
    使用tabix对vcf文件进行bgzip压缩
    """
    runshell(f"bgzip -f {vcf_file}")
    runshell(f"tabix -f -p vcf {vcf_file}.gz")


## outside command 
def replaceVCFsampleID(args):
    """
    Replace sample ID in VCF file
    >>> %(prog)s -s <replace_file.txt> <VCF.gz> [Options]
    -s <replace_file.txt> is formatted as follows:
    -----------------------------------------------
    old_id1 new_id1
    old_id2 new_id2
    ...
    -----------------------------------------------
    Note:The replace_file has no header by default!
    """ 
    install()
    p = argparse.ArgumentParser(prog=replaceVCFsampleID.__name__,
                        description=replaceVCFsampleID.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-s', '--replace_file', required=True,
                        help='Input the replace file')
    pReq.add_argument('vcf_file',
                        help='Input the replace vcf file')
    pReq.add_argument('-o','--output',
            help='Output the replaced VCF file')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    check_file_exists(args.replace_file)
    check_file_exists(args.vcf_file)

    if args.output.endswith('.gz'):
        output = args.output.replace('.gz','')
        log.info(f'Recognizes that the output is specified as a .gz file, and the output before compression is specified as [green blink]{output}[/]', extra={"markup": True})
    else:
        output = args.output

    replace_dict = read_rename_file(args.replace_file)
    rename_vcf_header(args.vcf_file,replace_dict,output)
    if os.path.exists(output):
        bgzip_vcf(output)
        log.info(f"Replace sample ID in VCF file [bold green blink]{output}[/] successfully", extra={"markup": True})
        log.info(f"Bgzip file: {output} to [bold green blink]{output}.gz[/]", extra={"markup": True})
    else:
        log.error(f"Replace sample ID in VCF file [bold red blink]{output}[/] failed!", extra={"markup": True})


def main():
    actions = (
            ("replaceVCFsampleID", "Replace sample ID in VCF file"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == '__main__':
    main()