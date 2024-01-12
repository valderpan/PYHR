#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/05/31

import re
import sys
import argparse
from rich import print
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify, read_file


log = richlog()

def parse_bamcov(covfile):
    Chr2cov = {}
    with open(covfile) as f:
        lines = (line.strip() for line in f)
        for line in lines:
            line_list = line.split()
            if int(line_list[3]) != 0:
                if line_list[0] not in Chr2cov:
                    Chr2cov[line_list[0]] = 0
                    Chr2cov[line_list[0]] += int(line_list[2]) - int(line_list[1])
                else:
                    Chr2cov[line_list[0]] += int(line_list[2]) - int(line_list[1])
    return Chr2cov


def read_fai(file):
    df = read_file(file,['seqid','len','uk1','uk2','uk3'])
    Chr2len = df.set_index('seqid').to_dict()['len']
    return Chr2len


def Cov_count(Chr2cov,Chr2len,outputfile=None):
    Gen_len = 0 ;Cov_len = 0
    print('{}\t{}\t{}\t{}'.format('seqid','Cov_len','Gen_cov','ratio'),file=outputfile)
    for key in Chr2cov.keys():
        print('{}\t{}\t{}\t{}'.format(key,Chr2cov[key],Chr2len[key],round(Chr2cov[key]/Chr2len[key],5)*100),file=outputfile)
        Gen_len += Chr2len[key]
        Cov_len += Chr2cov[key]
    print('{}\t{}\t{}\t{}'.format('Total',Cov_len,Gen_len,round(Cov_len/Gen_len,5)*100),file=outputfile)


def depth_count(depthfile,outputfile=None): 
    reads_num = 0
    coverage_depth = 0
    contig = ''
    with open(depthfile) as f:
        for line in f:
            line_list = line.split('\t')
            if reads_num == 0:
                coverage_depth += int(line_list[2].rstrip())
                reads_num += 1
                contig = line_list[0].rstrip()
            elif reads_num != 0 and line_list[0] == contig:
                coverage_depth += int(line_list[2].rstrip())
                reads_num += 1
            else:
                depth = float(coverage_depth) / float(reads_num)
                reads_num = 0
                coverage_depth = 0
                output = contig + '\t' + str(depth)
                print(output.rstrip(),file=outputfile)
        last_depth = float(coverage_depth)/float(reads_num)
        print(contig+'\t'+str(last_depth),file=outputfile)


def read_result(files):
    file2rowL = {}
    for file in files:
        with open(file) as f:
            rowL = [line.strip() for line in f]
        file2rowL[file] = rowL
    return file2rowL


def re_findends(patten,rowL):
    for i in rowL:
        if i.endswith(patten):
            return int(re.findall('([0-9]+)',i)[0])
        else:
            continue


def parse_result(rowL):
    single_reads = re_findends('reads; of these:',rowL)
    uniq_reads1 = re_findends('aligned concordantly exactly 1 time',rowL)
    uniq_reads2 = re_findends('aligned discordantly 1 time',rowL)
    uniq_reads3 = re_findends('aligned exactly 1 time',rowL)
    multi_reads1 = re_findends('aligned concordantly >1 times',rowL)
    multi_reads2 = re_findends('aligned >1 times',rowL)
    return [single_reads,uniq_reads1,uniq_reads2,uniq_reads3,multi_reads1,multi_reads2]


def cal_rate(reads_num_list):
    Total_reads = reads_num_list[0]*2
    uniq_reads = reads_num_list[1]*2+reads_num_list[2]*2+reads_num_list[3]
    uniq_rate = uniq_reads/Total_reads*100
    multi_reads = reads_num_list[4]*2+reads_num_list[5]
    multi_rate = multi_reads/Total_reads*100
    return [Total_reads,uniq_reads,uniq_rate,multi_reads,multi_rate]


def print_output(logfile,resL,output=None):
    print('Sample file : [bold magenta]{}[/bold magenta]'.format(logfile),file=output)
    print('    Total reads num : {}'.format(resL[0]),file=output)
    print('    Uniq Mapped reads num : {}'.format(resL[1]),file=output)
    print('    Uniq Mapped reads rate : {}%'.format(round(resL[2],2)),file=output)
    print('    Multi Mapped reads num : {}'.format(resL[3]),file=output)
    print('    Multi Mapped reads rate : {}%'.format(round(resL[4],2)),file=output)
    print('    Total Mapped reads num : {}'.format(resL[1]+resL[3]),file=output)
    print('    Total Mapped reads rate : {}%'.format(round((resL[1]+resL[3])/resL[0],4)*100),file=output)


#outside command
def StatSeqCoverage(args):
    '''
    Count the sequencing coverage on each chromosome according to the results of `bedtools genomecov`
    >>> bedtools genomecov -bga -pc -ibam ${bam} > ${bam}.cov
    >>> %prog ${bam}.cov [-o ${bam}.cov.stat] [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=StatSeqCoverage.__name__,
                        description=StatSeqCoverage.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')


    pReq.add_argument('coverage',
            help='Input the result file generated by `bedtools genomecov`')
    pReq.add_argument('fai',
            help='Input the result file generated by `samtools faidx`')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            #default=sys.stdout, 
            help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.coverage)
    check_file_exists(args.fai)
    Chr2COV = parse_bamcov(args.coverage)
    Chr2Len = read_fai(args.fai)

    if args.output:
        Cov_count(Chr2COV,Chr2Len,args.output)
        log.info('Completed! the output file is `{}`'.format(args.output.name))
    else:
        Cov_count(Chr2COV,Chr2Len)


def StatSeqDepth(args):
    '''
    Counting the sequencing depth of each chromosome according to the results of `samtools depth`
    >>> samtools depth -aa bam > depth.info
    >>> %prog depth.info [-o depth.info.stat] [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=StatSeqDepth.__name__,
                        description=StatSeqDepth.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')


    pReq.add_argument('depth',
            help='Input the result file generated by `samtools depth`')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            #default=sys.stdout, 
            help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.depth)

    if args.output :
        depth_count(args.depth,args.output)
        log.info('Completed! the output file is `{}`'.format(args.output.name))
    else:
        depth_count(args.depth)


def StatMappedRate(args):
    '''
    Calculate the comparison of each indicator after hisat2 runs by specifying the results obtained with the --summary-file parameter or a file with output directed by `&>`(such as bowtie2)
    Presently Known Supported Software : [hisat2]/[bowtir2]
    >>> %(prog)s <stat.file> [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=StatMappedRate.__name__,
                        description=StatMappedRate.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')


    pReq.add_argument('stat',nargs='+',
            help='Input the statistics file specified by hisat2')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    for file in args.stat:
        check_file_exists(file)
    file2rowL = read_result(args.stat)

    if args.output :
        for file in file2rowL:
            readsL = parse_result(file2rowL[file])
            resL = cal_rate(readsL)
            print_output(file,resL,args.output)
            log.info('Completed! The output file is `{}` '.format(args.output.name))
    else:
        for file in file2rowL:
            readsL = parse_result(file2rowL[file])
            resL = cal_rate(readsL)
            print_output(file,resL)


def main():
    actions = (
            ("StatSeqCoverage", "Count the sequencing coverage on each chromosome according to the results of `bedtools genomecov`"),
            ("StatSeqDepth", "Counting the sequencing depth of each chromosome according to the results of `samtools depth`"),
            ("StatMappedRate", "Calculate the comparison of each indicator after hisat2 runs by specifying the results obtained with the --summary-file parameter"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()


