#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/09/08


import os
import re
import sys
import argparse
import subprocess
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

class md5():
    def __init__(self,file) -> None:
        self.md5D = {}
        self.file = file
    
    def set_md5D(self):
        with open(self.file) as f:
            lines = (line.strip() for line in f)
            for line in lines:
                md5_value,file_name = line.split()
                if '/' in file_name:
                    file_name = file_name.split('/')[-1]
                self.md5D[file_name] = md5_value
        return self.md5D


def compare_md5(Omd5D,Nmd5D,pattern):
    all_files = len(Omd5D.keys())
    Error_match = 0
    Correct_match = 0
    if pattern == 'all':
        if len(Omd5D) != len(Nmd5D):
            log.error('The number of files before and after the transfer is different, please check it !!!')
            sys.exit()
        else:
            log.debug('Consistent number of files before and after transfer ~')
            for key in Omd5D.keys():
                if not key in Nmd5D.keys():
                    # log.error('{} md5 value is missing'.format(key))
                    log.warning('{} md5 value is missing'.format(key))
                else:
                    if Nmd5D[key] == Omd5D[key]:
                        log.info('{} md5 value is OK ~'.format(key))
                        Correct_match += 1
                    else:
                        log.error('[underline yellow]{} [/]md5 value is incorrectly checked !!!'.format(key),extra={"markup": True})
                        log.debug('The original md5 value of file {} is [underline green]{}[/]'.format(key,Omd5D[key]),extra={"markup": True})
                        log.warning('The transferred md5 value of file {} is [underline yellow]{}[/]'.format(key,Nmd5D[key]),extra={"markup": True})
                        Error_match += 1
        log.info('#--------------------------------------------------------------------#')
        log.info(f'Total transferred file number : {len(Nmd5D.keys())}')
        log.info('Correct checked file number : {}, Incorrect checked file number is {}'.format(Correct_match,Error_match))
        log.info('#=========================MD5 SUMMARY================================#')
        if Correct_match == all_files:
            log.info('#                 [bold green blink] All files are correctly checked ~ [/]                #', extra={"markup": True})
        else:
            log.warning('#     [bold red] Not all files are correctly checked, please check it !!! [/]     #', extra={"markup": True})
        log.info('#====================================================================#')
    elif pattern == 'part':
        if len(Omd5D) != len(Nmd5D):
            log.warning('The number of files before and after the transfer is different, but [-p=part] is triggered ,so it will still run!')
        else:
            log.warning('Consistent number of files before and after transfer, please specify the parameter pattern as [bold magenta blink]all[/] pattern [underline magenta blink](just use -p all)[/] ',extra={"markup": True})
        for key in Omd5D.keys():
            if not key in Nmd5D.keys():
                log.warning('`{}` md5 value is missing'.format(key))
            else:
                if Nmd5D[key] == Omd5D[key]:
                    log.info('`{}` md5 value is OK ~'.format(key))
                    Correct_match += 1
                else:
                    log.error('`{}` md5 value is incorrectly checked'.format(key))
                    log.debug('The original md5 value of file {} is {}'.format(key,Omd5D[key]))
                    log.warning('The transferred md5 value of file {} is {}'.format(key,Nmd5D[key]))
                    Error_match += 1
        log.info('#--------------------------------------------------------------------#')
        log.info(f'Total original file number : {len(Omd5D.keys())}')
        log.info(f"Total transferred file number : {len(Nmd5D.keys())}")
        log.info('Correct checked file number : {}, Incorrect checked file number is {}'.format(Correct_match,Error_match))
        log.info('#=========================MD5 SUMMARY================================#')
        if Correct_match == len(Nmd5D.keys()):
            log.info('#                 [bold green blink] All files are correctly checked ~ [/]                #', extra={"markup": True})
        else:
            log.warning('#     [bold red] Not all files are correctly checked, please check it !!! [/]     #', extra={"markup": True})
        log.info('#====================================================================#')


def download_fastq(SRR_file):
    total_file = 0
    with open(SRR_file) as f:
        lines = (line.strip() for line in f)
        for line in lines:
            if not line.startswith('SRR'):
                continue
            else:
                total_file += 1
                line_list = line.split()
                if os.path.exists(f"{line_list[1]}.sra") or os.path.exists(f"{line_list[1]}.fastq.gz"):
                    log.info(f"[bold green blink] {line_list[1]}.sra/fastq.gz already exists, skip![/]",extra={"markup": True})
                else:
                    log.info(f'Download [bold green blink]{line_list[0]} [/]and store it in file [bold green blink] {line_list[1]}.sra [/]',extra={"markup": True})
                    cmd1 = f'prefetch --max-size 200G {line_list[0]} -o {line_list[1]}.sra'
                    runshell(cmd1)
                    cmd2 = f'pfastq-dump -t 10 --gzip --split-3 {line_list[1]}.sra'
                    log.info(f'Convert {line_list[1]}.sra to {line_list[1]}.fastq.gz')
                    runshell(cmd2)

    files_num = len([file for file in Path('./').files() if file.endswith('sra')])
    if total_file != files_num:
        log.error('Not all files are downloaded, please check it!')


def download_GSE_suppl(accession,mode,filename = None,output_dir = './'):
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        log.info(f"创建输出目录: {output_dir}")
    
    if mode == 'GSE' and accession.startswith('GSE'):
        #获取GSE accession的前6位,GSE加上前三位数字
        prefix6 = accession[:6]
        cmd1 = f"wget -P {output_dir} ftp://ftp.ncbi.nlm.nih.gov/geo/series/{prefix6}nnn/{accession}/suppl/{accession}_RAW.tar"
        log.info(f"正在下载GSE数据: {accession}")
        runshell(cmd1)
        output_path = os.path.join(output_dir, f"{accession}_RAW.tar")
        # 修改4: 添加log.info打印下载文件存储位置
        log.info(f"GSE数据已下载到: {output_path}")

    elif mode == 'GSM' and accession.startswith('GSM'):
        #下载单个文件
        #首先判断传入的accession是否是GSM accession,如果是GSM accession则直接下载对应的文件
        if filename is None:
            log.error(f"GSM模式需要指定filename参数")
            return False
        
        prefix7 = accession[:7]
        cmd2 = f"wget -P {output_dir} ftp://ftp.ncbi.nlm.nih.gov/geo/samples/{prefix7}nnn/{accession}/suppl/{filename}"
        log.info(f"正在下载GSM数据: {accession}/{filename}")
        runshell(cmd2)
        output_path = os.path.join(output_dir, filename)
        # 修改7: 添加log.info打印下载文件存储位置
        log.info(f"GSM数据已下载到: {output_path}")

    else:
        log.error('Option {} does not exist, it can only be one of [GSE or GSM]'.format(mode))



def mv_encode_seq(metadata,seqpath):
    df = pd.read_table(metadata,sep='\t')
    biosample = df['Biosample term name'].unique()
    for s in biosample:
        qdf = df[df['Biosample term name'] == s]
        ENC_D = {}
        for index,row in qdf.iterrows():
            paired = row['Paired with'].split('/')[2]
            if row['File accession'] not in ENC_D and row['File accession'] not in list(ENC_D.values()):
                ENC_D[row['File accession']] = paired
        for num,key in enumerate(ENC_D.keys()):
            cmd1 = f'mv {seqpath}/{key}.fastq.gz {seqpath}/{s}.Lib{num+1}.R1.fastq.gz'
            cmd2 = f'mv {seqpath}/{ENC_D[key]}.fastq.gz {seqpath}/{s}.Lib{num+1}.R2.fastq.gz'
            runshell(cmd1)
            runshell(cmd2)


def EvaluatingSeqDepth(seqpath,output):
    Seqfiles = [file for file in Path(seqpath).files() if file.endswith('.fq.gz') or file.endswith('.fastq.gz')]
    file2num = {}
    pairedD = {}
    for file in Seqfiles:
        cmd1 = f'zcat {file} | echo $((`wc -l`/4))'
        log.info(f'Run commond {cmd1}')
        res = subprocess.check_output(cmd1, shell=True)
        res = res.decode('utf-8').strip()
        file2num[str(file.basename())] = res
    for key,value in file2num.items():
        if value in pairedD:
            pairedD[value].append(key)
        else:
            pairedD[value] = [key]
    for key,value in pairedD.items():
        if len(value) > 1:
            base = int(key)*2*150
            depth = round(base/(3*(10**9)),1)
            # print('+'.join(value),base,depth,sep='\t',file=output)
        else:
            base = int(key)*150
            depth = round(base/(3*(10**9)),1)
            # print('+'.join(value),base,depth,sep='\t',file=output)
    table = Table(show_header=True, header_style="bold magenta")
    table.add_column("File Name", style="cyan", width=50)
    table.add_column("Base Num",style='green')
    table.add_column("SequencingDepth", justify="right",style='magenta')
    with open(output,'w') as w1:
        for key,value in pairedD.items():
            if len(value) > 1:
                base = int(key)*2*150
                depth = round(base/(3*(10**9)),1)
                # print('+'.join(value),base,depth,sep='\t',file=output)
                w1.write(f"{'+'.join(value)}\t{base}\t{depth}\n")
                table.add_row('+'.join(value),str(base),str(depth))
            else:
                base = int(key)*150
                depth = round(base/(3*(10**9)),1)
                # print('+'.join(value),base,depth,sep='\t',file=output)
                w1.write(f"{'+'.join(value)}\t{base}\t{depth}\n")
                table.add_row('+'.join(value),str(base),str(depth))
    console.print(table)


## outside command 
def CheckMD5(args):
    """
    Correct the md5 value before and after file transfer
    >>> %(prog)s <Original md5 file> <transferred md5 file> [Options]
    -p [pattern] means Whether the new md5 file after transmission is part or all of the Original md5 file
    """ 
    install()
    p = argparse.ArgumentParser(prog=CheckMD5.__name__,
                        description=CheckMD5.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('OriginalMD5',
            help='Input original md5 file')
    pReq.add_argument('TransferredMD5',
            help='Input the transferred file')
    pReq.add_argument('-p','--pattern',choices=['all','part'],
                      required=True,
            help='Input the transferred file')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    check_file_exists(args.OriginalMD5)
    check_file_exists(args.TransferredMD5)
    od = md5(args.OriginalMD5)
    nd = md5(args.TransferredMD5)
    pattern = args.pattern
    compare_md5(od.set_md5D(),nd.set_md5D(),pattern)


def DownloadFastq(args):
    """
    Download GEO database public data and convert it to .fastq.gz through python
    Attention: The SRR file must contain two columns:[SRRnumber sample_name],separated by space or tab !    
    >>> %(prog)s <SRR_with_name file> [Options]
    """ 
    install()
    p = argparse.ArgumentParser(prog=DownloadFastq.__name__,
                        description=DownloadFastq.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('SRRfile',
            help='Input SRR_with_name file')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    args = p.parse_args(args)

    check_file_exists(args.SRRfile)
    download_fastq(args.SRRfile)


def MVENCODE(args):
    """
    Rename the ENCODE data according to metadata.tsv
    Attention: The metadata.tsv file must contain three columns:[File accession,Paired with,Biosample term name]
    >>> %(prog)s <metadata.tsv> [seqpath] [Options]
    """ 
    install()
    p = argparse.ArgumentParser(prog=MVENCODE.__name__,
                        description=MVENCODE.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('metadata',
            help='Input metadata.tsv file')
    pReq.add_argument('seqpath',
            help='Input fastq path')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    args = p.parse_args(args)

    check_file_exists(args.metadata)
    mv_encode_seq(args.metadata,args.seqpath)


def EvaluateSeqDepth(args):
    """
    Evaluate the sequencing volume of the sequencing file
    Attention: Fastq sequencing files must end in.fq.gz or.fastq.gz!!!
    >>> %(prog)s seqpath [-o output]
    """ 
    install()
    p = argparse.ArgumentParser(prog=EvaluateSeqDepth.__name__,
                        description=EvaluateSeqDepth.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('seqpath',
            help='Input fastq path')
    pReq.add_argument('-o', '--output', required=True,
                        help='output file ')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    EvaluatingSeqDepth(args.seqpath,args.output)


## outside command
def DownloadGSESuppl(args):
    """
    Download the supplementary files (RAW.tar or single files) from the GEO database.

    Examples:
    1. Download a GSE's RAW.tar file:
        %(prog)s GSE12345 -m GSE -o ./downloads
    2. Download a specific supplementary file for a GSM:
        %(prog)s GSM1234567 -m GSM -f GSM1234567_file.txt -o ./gsm_files

    >>> %(prog)s <accession> <mode> [Options]
    """
    install()
    p = argparse.ArgumentParser(prog=DownloadGSESuppl.__name__,
                        description=DownloadGSESuppl.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('accession',
                       help='GEO accession number (e.g., GSE12345 or GSM1234567)')
    pReq.add_argument('-m','--mode', choices=['GSE', 'GSM'],
                       help='Download mode. "GSE" for a series, "GSM" for a sample.')
    pOpt.add_argument('-f', '--filename', default=None,
                       help='Filename to download (REQUIRED if mode is "GSM").')
    pOpt.add_argument('-o', '--output_dir', default='./',
                       help='Output directory for downloaded files (default: current directory).')
    pOpt.add_argument('-h', '--help', action='help',
                       help='show help message and exit.')

    args = p.parse_args(args)

    # 验证参数逻辑
    if args.mode == 'GSM' and not args.filename:
        p.error('The following argument is required when mode is GSM: -f/--filename')

    # 调用核心下载函数
    download_GSE_suppl(args.accession, args.mode, args.filename, args.output_dir)
    log.info(f"[bold green blink]Download completed successfully![/]", extra={"markup": True})


def main():
    actions = (
            ("CheckMD5", "Correct the md5 value before and after file transfer"),
            ("DownloadFastq","Download GEO database public data and convert it to .fastq.gz through python"),
            ("MVENCODE","Rename the ENCODE data according to metadata.tsv"),
            ("EvaluateSeqDepth","Evaluate the sequencing volume of the sequencing file"),
            ("DownloadGSESuppl","Download the supplementary files of GSE from GEO database"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()