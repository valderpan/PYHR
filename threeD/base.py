#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2024/03/13


import cooler
import bioframe
import cooltools
import numpy as np
import pandas as pd
import cooltools.lib.plotting
import matplotlib.pyplot as plt
from path import Path
from collections import OrderedDict
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify, read_file, runshell


log = richlog()

def read_cool(coolfile):
    """
    Input the balanced cool file
    The cool file can be obtained by converting HiC-Pro's raw matrix with hicConvertFormat, 
    and then normalizing the cool file with cooler balance
    """
    clr = cooler.Cooler(coolfile)
    return clr


def readGCfile(gcfile): 
    """
    Input the genome gc file
    In humans and mice, GC content is useful for phasing \
    because it typically has a strong correlation at the 100kb-1Mb bin level with the eigenvector. 
    In other organisms, other phasing tracks have been used to orient eigenvectors from Hi-C data.
    """
    gc_cov = pd.read_table(gcfile)
    return gc_cov


def getResolution(clr):
    """
    Get the resolution of the cool file
    """
    resolution = clr.info['bin-size']
    return resolution


def countGC(fasta,clr,output,output_Dir = '.'):
    """
    Count the GC content of genome
    """
    bins = clr.bins()[:]
    genome = bioframe.load_fasta(fasta);
    ## note the next command may require installing pysam
    gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], genome)
    gc_cov.to_csv(f'{output_Dir}/{output}',index=False,sep='\t')
    log.info(f'Completed! the GCcontent file is {output_Dir}/{output}')


def filter_chrom(clr):
    view_df = pd.DataFrame({'chrom': clr.chromnames,
                        'start': 0,
                        'end': clr.chromsizes.values,
                        'name': clr.chromnames}
                      )
    view_df_clean = view_df[~view_df['chrom'].str.contains('chrM')]
    view_df_clean = view_df_clean[view_df_clean['chrom'].str.contains('chr')]
    view_df_clean = view_df_clean[view_df_clean['chrom'].str.len() < 6]
    return view_df_clean


def calExpectedCis(clr,
    view_df,
    intra_only=True,
    smooth=True,
    aggregate_smoothed=True,
    smooth_sigma=0.1,
    clr_weight_name="weight",
    ignore_diags=2,  # should default to cooler info
    chunksize=10_000_000,
    nproc=1):
    """
    Copy from cooltools:https://github.com/open2c/cooltools/blob/master/cooltools/cli/expected_cis.py

    Calculate average interaction frequencies as a function of genomic
    separation between pixels i.e. interaction decay with distance.
    Genomic separation aka "dist" is measured in the number of bins,
    and defined as an index of a diagonal on which pixels reside (bin1_id - bin2_id).

    Average values are reported in the columns with names {}.avg, and they
    are calculated as a ratio between a corresponding sum {}.sum and the
    total number of "valid" pixels on the diagonal "n_valid".

    When balancing weights (clr_weight_name=None) are not applied to the data, there is no
    masking of bad bins performed.


    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object
    view_df : viewframe
        a collection of genomic intervals where expected is calculated
        otherwise expected is calculated for full chromosomes.
        view_df has to be sorted, when inter-regions expected is requested,
        i.e. intra_only is False.
    intra_only: bool
        Return expected only for symmetric intra-regions defined by view_df,
        i.e. chromosomes, chromosomal-arms, intra-domains, etc.
        When False returns expected both for symmetric intra-regions and
        assymetric inter-regions.
    smooth: bool
        Apply smoothing to cis-expected. Will be stored in an additional column
    aggregate_smoothed: bool
        When smoothing, average over all regions, ignored without smoothing.
    smooth_sigma: float
        Control smoothing with the standard deviation of the smoothing Gaussian kernel.
        Ignored without smoothing.
    clr_weight_name : str or None
        Name of balancing weight column from the cooler to use.
        Use raw unbalanced data, when None.
    ignore_diags : int, optional
        Number of intial diagonals to exclude results
    chunksize : int, optional
        Size of pixel table chunks to process
    nproc : int, optional
        How many processes to use for calculation. Ignored if map_functor is passed.
    map_functor : callable, optional
        Map function to dispatch the matrix chunks to workers.
        If left unspecified, pool_decorator applies the following defaults: if nproc>1 this defaults to multiprocess.Pool;
        If nproc=1 this defaults the builtin map. 

    Returns
    -------
    DataFrame with summary statistic of every diagonal of every symmetric
    or asymmetric block:
    region1, region2, diag, n_valid, count.sum count.avg, etc

    """

    # cvd == contacts-vs-distance
    cvd = cooltools.expected_cis(
            clr,
            view_df,
            intra_only=intra_only,
            smooth=smooth,
            aggregate_smoothed=aggregate_smoothed,
            smooth_sigma=smooth_sigma,
            clr_weight_name=clr_weight_name,
            ignore_diags=ignore_diags,  # should default to cooler info
            chunksize=chunksize,
            nproc=nproc)
    return cvd

