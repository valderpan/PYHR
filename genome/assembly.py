#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/05/16



#TODO  这个代码还没有正式开始写


import sys
import gzip
import time
import argparse
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.utils.fasta import Fa2dict


log = richlog()

def depth_count(depthfile,outputfile):
    '''
    计算每个contig的覆盖深度
    :param depthfile: samtools depth -aa -b contig1.scale *.sorted.bam > depthfile
    :param outputfile:
    :return:
    '''
    reads_num = 0
    coverage_depth = 0
    contig = ''
    result = open(outputfile,'w')
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
                result.write(output.rstrip()+'\n')
        last_depth = float(coverage_depth)/float(reads_num)
        result.write(contig+'\t'+str(last_depth)+'\n')
    result.close()
