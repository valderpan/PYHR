#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2023/12/04


import re
import sys
import cooler
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from coolpuppy import coolpup
from coolpuppy.lib.numutils import get_domain_score
from coolpuppy.lib.puputils import accumulate_values
from coolpuppy import plotpup
from cooltools.lib import io
from cooltools.lib import plotting
from path import Path
from natsort import natsorted
from collections import OrderedDict
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify, read_file
from PYHR.utils.fasta import Fa2dict
from PYHR.threeD.base import read_cool, filter_chrom, calExpectedCis

log = richlog()


def parseBoundaries(Boundaries_bed):
    with open(Boundaries_bed,'r') as f:
        prefix = re.findall('([0-9A-Za-z\.\_]+)\-\-',Boundaries_bed)[0]
        with open(f'{prefix}.TAD.bed','w') as w:
            lines = f.readlines()
            for i in range(0,len(lines)):
                if i == 0 or i+1 >= len(lines):
                    continue
                else:
                    fields1 = lines[i].strip().split('\t')
                    fields2 = lines[i+1].strip().split('\t')
                    w.write(f"{fields1[0]}\t{fields1[2]}\t{fields2[1]}\n")
    return prefix


def readTAD(TADfile):
    """
    TAD combined results based on Cword-matrix2insulation.pl
    The file contains only three columns,:[seqnames, bin_start, bin_end]
    """
    TAD = pd.read_table(TADfile, header=None, names=['chrom', 'start', 'end'])
    return TAD


def add_domain_score(snippet):
    """
    Define a helper function to store domain scores within each snippet
    """
    snippet['domain_score'] = get_domain_score(snippet['data']) # Calculates domain score for each snippet according to Flyamer et al., 2017
    return snippet


def extra_sum_func(dict1, dict2):
    """
    Another helper function to save domain scores when combining snippets into a pileup
    """
    return accumulate_values(dict1, dict2, 'domain_score')


def runATA(clr, tad, res, view_df, cvd):
    cc = coolpup.CoordCreator(tad, resolution=res, features_format='bed', local=True, rescale_flank=1)
    pu = coolpup.PileUpper(clr, cc, expected=cvd, view_df=view_df, ignore_diags=0, rescale_size=99, rescale=True)
    pup = pu.pileupsWithControl(postprocess_func=add_domain_score, # Any function can be applied to each snippet before they are averaged in the postprocess_func
                            extra_sum_funcs={'domain_score': extra_sum_func}) # If additional values produced by postprocess_func need to be saved,
                                                                              # it can be done using the extra_sum_funcs dictionary, which defines how to combine them.
    return pup


def outputATA(pup, prefix,output_Dir = '.'):
    pd.DataFrame(pup.loc[0,"domain_score"],columns=['TAD_score']).to_csv(f'{output_Dir}/ATA_TADscore.{prefix}.csv',index=False)
    pd.DataFrame(pup.loc[0,'data']).to_csv(f'{output_Dir}/ATA_rescale_matrix.{prefix}.csv',index=False)
    log.info(f'Completed! the ATA rescale matrix is {output_Dir}/ATA_rescale_matrix.{prefix}.csv')
    log.info(f'Completed! the ATA TADscore is {output_Dir}/ATA_TADscore.{prefix}.csv')


def plotATA(pup, prefix,vmin = None, vmax = None,output_Dir = '.'):
    plotpup.plot(pup,
             score=False,
             height=5,
             cmap = 'coolwarm',
             font='Arial',
             vmin=vmin,
             vmax=vmax)
    plt.savefig(f'{output_Dir}/ATA.plot.{prefix}.pdf', dpi=300, bbox_inches='tight')
    log.info(f'Completed! the ATA plot is {output_Dir}/ATA.plot.{prefix}.pdf')


## outside command
def boundaries2bed(args):
    """
    Convert boundaries to bed based on cworld results
    The result of running this script is equivalent to [Cworld-insulation2tads.pl] !!!
    >>> %(prog)s <boundaries.bed> [Options]
    """ 
    install()
    p = argparse.ArgumentParser(prog=boundaries2bed.__name__,
                        description=boundaries2bed.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    
    pReq.add_argument('BoundariesBed',
                    type=str,
                    help="Input cworld's result file boundaries.bed")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.BoundariesBed)
    res_prefix = parseBoundaries(args.BoundariesBed)
    log.info(f'The results file is output to `{res_prefix}.TAD.bed`')


def ATA(args):
    """
    ATA analysis based on coolpuppy
    >>> %(prog)s -c <cool file> -t <tad file> -n <output name> [Options]
    """
    p = argparse.ArgumentParser(prog=ATA.__name__,
                        description=ATA.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('-c','--cool',
                    type=str,
                    help="Input the balanced cool file")
    pReq.add_argument('-t','--tad',
                    type=str,
                    help="Input the TAD file")
    pReq.add_argument('-n','--name',
                    type=str,
                    help="Input the prefix of the output file")
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    pOpt.add_argument('-vmin',
            help='Specifies the minimum value of the legend in the ATA plot')
    pOpt.add_argument('-vmax',
            help='Specifies the maximum value of the legend in the ATA chart')
    pOpt.add_argument('-d','--output_Dir',default='.',
                      help='The output directory of the output file')

    args = p.parse_args(args)
    check_file_exists(args.cool)
    check_file_exists(args.tad)
    clr = read_cool(args.cool)
    tad = readTAD(args.tad)
    view_df = filter_chrom(clr)
    cvd = calExpectedCis(clr,ignore_diags=0,view_df = view_df,chunksize=1000000)
    pup = runATA(clr, tad, res=clr.info['bin-size'], view_df=view_df, cvd=cvd)
    outputATA(pup, args.name, args.output_Dir)
    plotATA(pup, args.name, args.vmin, args.vmax, args.output_Dir)


def main():
    actions = (
            ("boundaries2bed", "Convert boundaries to bed based on cworld results"),
            ("ATA", "ATA analysis based on coolpuppy"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()