#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2023/01/21

import re
import sys
import math
import cooler
import bioframe
import argparse
import cooltools
import numpy as np
import pandas as pd
import cooltools.lib.plotting
import matplotlib.pyplot as plt
from path import Path
from scipy import sparse
from cytoolz import merge
from collections import OrderedDict
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify, read_file, runshell


log = richlog()
#Set the parameters for saddle plot
Q_LO = 0.025 # ignore 2.5% of genomic bins with the lowest E1 values
Q_HI = 0.975 # ignore 2.5% of genomic bins with the highest E1 values
N_GROUPS = 38 # divide remaining 95% of the genome into 38 equisized groups, 2.5% each


def read_cool(coolfile):
    """
    Input the balanced cool file
    The cool file can be obtained by converting HiC-Pro's raw matrix with hicConvertFormat, 
    and then normalizing the cool file with cooler balance
    """
    clr = cooler.Cooler(coolfile)
    return clr


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


def readGCfile(gcfile): 
    """
    Input the genome gc file
    In humans and mice, GC content is useful for phasing \
    because it typically has a strong correlation at the 100kb-1Mb bin level with the eigenvector. 
    In other organisms, other phasing tracks have been used to orient eigenvectors from Hi-C data.
    """
    gc_cov = pd.read_table(gcfile)
    return gc_cov


def filter_chrom(clr):
    view_df = pd.DataFrame({'chrom': clr.chromnames,
                        'start': 0,
                        'end': clr.chromsizes.values,
                        'name': clr.chromnames}
                      )
    view_df_clean = view_df[~view_df['chrom'].str.contains('chrM')]
    view_df_clean = view_df_clean[view_df_clean['chrom'].str.contains('chr')]
    return view_df_clean


def callCisEigs(clr,gc_cov,view_df):
    # obtain first 3 eigenvectors
    cis_eigs = cooltools.eigs_cis(
                            clr,
                            gc_cov,
                            view_df=view_df,
                            n_eigs=3,
                            )
    # cis_eigs[0] returns eigenvalues, here we focus on eigenvectors
    eigenvector_track = cis_eigs[1][['chrom','start','end','E1']]
    return eigenvector_track


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


def calSaddle(clr,cvd,eigenvector_track,view_df):
    interaction_sum, interaction_count =  cooltools.saddle(
        clr,
        cvd,
        eigenvector_track,
        'cis',
        n_bins=N_GROUPS,
        qrange=(Q_LO,Q_HI),
        view_df=view_df)
    return interaction_sum, interaction_count


def outputSum_Count(interaction_sum,interaction_count,prefix,output_Dir = '.'):
    pd.DataFrame(interaction_sum).to_csv(f'{output_Dir}/interaction_sum.{prefix}.csv',
                        index=False)
    pd.DataFrame(interaction_count).to_csv(f'{output_Dir}/interaction_count.{prefix}.csv',
                        index=False)
    log.info(f'Completed! the interaction file is {output_Dir}/interaction_sum.{prefix}.csv and {output_Dir}/interaction_count.{prefix}.csv')


def saddleplot(
    track,
    saddledata,
    n_bins,
    vrange=None,
    qrange=(0.0, 1.0),
    cmap="coolwarm",
    scale="log",
    vmin=0.5,
    vmax=2,
    color=None,
    title=None,
    xlabel=None,
    ylabel=None,
    clabel=None,
    fig=None,
    fig_kws=None,
    heatmap_kws=None,
    margin_kws=None,
    cbar_kws=None,
    subplot_spec=None):
    """
    Copy from cooltools::https://github.com/open2c/cooltools/blob/master/cooltools/api/saddle.py

    Generate a saddle plot.
    Parameters
    ----------
    track : pd.DataFrame
        See cooltools.digitize() for details.
    saddledata : 2D array-like
        Saddle matrix produced by `make_saddle`. It will include 2 flanking
        rows/columns for outlier signal values, thus the shape should be
        `(n+2, n+2)`.
    cmap : str or matplotlib colormap
        Colormap to use for plotting the saddle heatmap
    scale : str
        Color scaling to use for plotting the saddle heatmap: log or linear
    vmin, vmax : float
        Value limits for coloring the saddle heatmap
    color : matplotlib color value
        Face color for margin bar plots
    fig : matplotlib Figure, optional
        Specified figure to plot on. A new figure is created if none is
        provided.
    fig_kws : dict, optional
        Passed on to `plt.Figure()`
    heatmap_kws : dict, optional
        Passed on to `ax.imshow()`
    margin_kws : dict, optional
        Passed on to `ax.bar()` and `ax.barh()`
    cbar_kws : dict, optional
        Passed on to `plt.colorbar()`
    subplot_spec : GridSpec object
        Specify a subregion of a figure to using a GridSpec.
    Returns
    -------
    Dictionary of axes objects.
    """

#     warnings.warn(
#         "Generating a saddleplot will be deprecated in future versions, "
#         + "please see https://github.com/open2c_examples for examples on how to plot saddles.",
#         DeprecationWarning,
#     )

    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from matplotlib.colors import Normalize, LogNorm
    from matplotlib import ticker
    import matplotlib.pyplot as plt

    class MinOneMaxFormatter(ticker.LogFormatter):
        def set_locs(self, locs=None):
            self._sublabels = set([vmin % 10 * 10, vmax % 10, 1])

        def __call__(self, x, pos=None):
            if x not in [vmin, 1, vmax]:
                return ""
            else:
                return "{x:g}".format(x=x)

    track_value_col = track.columns[3]
    track_values = track[track_value_col].values

    digitized_track, binedges = cooltools.digitize(
        track, n_bins, vrange=vrange, qrange=qrange
    )
    x = digitized_track[digitized_track.columns[3]].values.astype(int).copy()
    x = x[(x > -1) & (x < len(binedges) + 1)]

    # Old version
    # hist = np.bincount(x, minlength=len(binedges) + 1)

    groupmean = track[track.columns[3]].groupby(digitized_track[digitized_track.columns[3]]).mean()

    if qrange is not None:
        lo, hi = qrange
        binedges = np.linspace(lo, hi, n_bins + 1)

    # Barplot of mean values and saddledata are flanked by outlier bins
    n = saddledata.shape[0]
    X, Y = np.meshgrid(binedges, binedges)
    C = saddledata
    if (n - n_bins) == 2:
        C = C[1:-1, 1:-1]
        groupmean = groupmean[1:-1]

    # Layout
    if subplot_spec is not None:
        GridSpec = partial(GridSpecFromSubplotSpec, subplot_spec=subplot_spec)
    grid = {}
    gs = GridSpec(
        nrows=3,
        ncols=3,
        width_ratios=[0.2, 1, 0.1],
        height_ratios=[0.2, 1, 0.1],
        wspace=0.05,
        hspace=0.05,
    )

    # Figure
    if fig is None:
        fig_kws_default = dict(figsize=(5, 5))
        fig_kws = merge(fig_kws_default, fig_kws if fig_kws is not None else {})
        fig = plt.figure(**fig_kws)

    # Heatmap
    if scale == "log":
        norm = LogNorm(vmin=vmin, vmax=vmax)
    elif scale == "linear":
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError("Only linear and log color scaling is supported")

    grid["ax_heatmap"] = ax = plt.subplot(gs[4])
    heatmap_kws_default = dict(cmap="coolwarm", rasterized=True)
    heatmap_kws = merge(
        heatmap_kws_default, heatmap_kws if heatmap_kws is not None else {}
    )
    img = ax.pcolormesh(X, Y, C, norm=norm, **heatmap_kws)
    plt.gca().yaxis.set_visible(False)

    # Margins
    margin_kws_default = dict(edgecolor="k", facecolor=color, linewidth=1)
    margin_kws = merge(margin_kws_default, margin_kws if margin_kws is not None else {})
    # left margin hist
    grid["ax_margin_y"] = plt.subplot(gs[3], sharey=grid["ax_heatmap"])

    plt.barh(
        binedges, height=1/len(binedges), width=groupmean, align="edge", **margin_kws
    )

    plt.xlim(plt.xlim()[1], plt.xlim()[0])  # fliplr
    plt.ylim(hi, lo)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["bottom"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    # top margin hist
    grid["ax_margin_x"] = plt.subplot(gs[1], sharex=grid["ax_heatmap"])

    plt.bar(
        binedges, width=1/len(binedges), height=groupmean, align="edge", **margin_kws
    )

    plt.xlim(lo, hi)
    # plt.ylim(plt.ylim())  # correct
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    plt.gca().yaxis.set_visible(False)

    # Colorbar
    grid["ax_cbar"] = plt.subplot(gs[5])
    cbar_kws_default = dict(fraction=0.8, label=clabel or "")
    cbar_kws = merge(cbar_kws_default, cbar_kws if cbar_kws is not None else {})
    if scale == "linear" and vmin is not None and vmax is not None:
        grid["ax_cbar"] = cb = plt.colorbar(img, **cbar_kws)
        # cb.set_ticks(np.arange(vmin, vmax + 0.001, 0.5))
        # # do linspace between vmin and vmax of 5 segments and trunc to 1 decimal:
        decimal = 10
        nsegments = 5
        cd_ticks = np.trunc(np.linspace(vmin, vmax, nsegments) * decimal) / decimal
        cb.set_ticks(cd_ticks)
    else:
        print('cbar')

        cb = plt.colorbar(img, format=MinOneMaxFormatter(), cax=grid["ax_cbar"], **cbar_kws)
        cb.ax.yaxis.set_minor_formatter(MinOneMaxFormatter())

    # extra settings
    grid["ax_heatmap"].set_xlim(lo, hi)
    grid["ax_heatmap"].set_ylim(hi, lo)
    grid['ax_heatmap'].grid(False)
    if title is not None:
        grid["ax_margin_x"].set_title(title)
    if xlabel is not None:
        grid["ax_heatmap"].set_xlabel(xlabel)
    if ylabel is not None:
        grid["ax_margin_y"].set_ylabel(ylabel)

    return grid


def outputEigenvector(eigenvector_track,prefix,output_Dir = '.'):
    eigenvector_track.to_csv(f'{output_Dir}/Eigenvector.{prefix}.bed',
                             header=False,index=False,sep='\t')
    log.info(f'Completed! the Eigenvector file is {output_Dir}/Eigenvector.{prefix}.bed')


def readInteractionSum(interaction_sum):
    return pd.read_csv(interaction_sum).values


def readInteractionCount(interaction_count):
    return pd.read_csv(interaction_count).values


def readEigenvector(eigenvector_track):
    return pd.read_csv(eigenvector_track,sep='\t',names=['chrom','start','end','E1'])


def plotsaddle(interaction_sum,interaction_count,eigenvector_track,prefix,output_Dir = '.'):
    saddleplot(eigenvector_track,
           interaction_sum/interaction_count,
           N_GROUPS,
           qrange=(Q_LO,Q_HI),
           cbar_kws={'label':f'{prefix} average observed/expected contact frequency'},
           vmin=0.5,vmax=2
          )
    #save plot
    plt.savefig(f'{output_Dir}/saddleplot_{prefix}.pdf', dpi=300, bbox_inches='tight')
    log.info(f'Completed! the saddle plot is {output_Dir}/saddleplot_{prefix}.pdf')


def calSaddleStrength(interaction_sum,interaction_count):
    from cooltools.api.saddle import saddle_strength
    # at extent=0, this reduces to ((S/C)[0,0] + (S/C)[-1,-1]) / (2*(S/C)[-1,0])
    strength = saddle_strength(interaction_sum, interaction_count)
    return(strength)


def outputSaddleStrength(strength,prefix,output_Dir = '.'):
    pd.DataFrame(strength).to_csv(
        f'{output_Dir}/saddle_strength.{prefix}.csv',index=False,header=False
    )
    log.info(f'Completed! the SaddleStrength file is {output_Dir}/saddle_strength.{prefix}.csv')


def store_fpkm(fpkmdf):
    fpkmD = {}
    col_names = fpkmdf.columns.tolist()
    for i in range(1,fpkmdf.shape[1]):
        fpkmD[col_names[i]] = fpkmdf.set_index(col_names[0]).to_dict()[col_names[i]]
    return fpkmD


def assign_fpkm(beddf,fpkmD,prefix):
    for key in fpkmD.keys():
        N = 0
        concatL = []
        for index,row in beddf.iterrows():
            if not isinstance(row['ID'],float):
                fpkmL = []
                idL = row['ID'].split(',')
                for id in idL:
                    if id not in fpkmD[key].keys():
                        N +=1
                    else:
                        fpkmL.append(math.log2(fpkmD[key][id]+1))
                if len(fpkmL) == 0:
                    bin_mean = 0
                else:
                    bin_mean = sum(fpkmL)/len(fpkmL)
                new_row = pd.DataFrame({'seqname':[row['seqname']],
                'start':[row['start']],'end':[row['end']],
                'value':[bin_mean]})
                concatL.append(new_row)
            else:
                new_row = pd.DataFrame({'seqname':[row['seqname']],
                'start':[row['start']],'end':[row['end']],
                'value':np.nan})
                concatL.append(new_row)
        res = pd.concat(concatL)
        res.to_csv(f'{prefix}.{key}.FPKM.csv',header=False,index=False)
        log.info(f'Completed! the output file is {prefix}.{key}.FPKM.csv')
        print(f'Total missing number of sample {key} is {N}')


#outside command 
def CalBinFPKM(args):
    '''
    Calculate the average expression value of each bin
    >>> %(prog)s -f <fasta> -w windows_file [-o output] [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=CalBinFPKM.__name__,
                        description=CalBinFPKM.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('-e','--expmatrix',required=True,
            help='Input fasta file')
    pReq.add_argument('-b','--bed',required=True,
            help='Input the bin file containing the gene ID (available from utils.dist SlidingWindow2GeneID)')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    check_file_exists(args.expmatrix)
    fpkmdf = read_file(args.expmatrix)
    beddf = read_file(args.bed,['seqname','start','end','ID'])
    prefix = re.findall('([\.\_0-9a-zA-Z]+)\.bed',args.bed)[0]
    fpkmD = store_fpkm(fpkmdf)
    assign_fpkm(beddf,fpkmD,prefix)


def CountsGC(args):
    """
    Counting GC content of each bin based on the bin of Cool and fasta file
    >>> %(prog)s -f <fasta> -c cool_file -o output [Options]
    """
    install()
    p = argparse.ArgumentParser(prog=CountsGC.__name__,
                        description=CountsGC.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('-f','--fasta',required=True,
                help='Input fasta file')
    pReq.add_argument('-c','--cool',required=True,
                help='Input cool file')
    pReq.add_argument('-o','--output',required=True,
                      help='The prefix of the output file(i.e. Mm10.gc_cov.tsv)')
    pOpt.add_argument('-d','--output_Dir',default='.',
                      help='The output directory of the output file')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    check_file_exists(args.fasta)
    check_file_exists(args.cool)
    clr = read_cool(args.cool)
    countGC(args.fasta,clr,args.output,args.output_Dir)


def cooltoolsAB(args):
    '''
    Call compartment based on cooltools algorithm
    >>> %(prog)s -c <cool_file> -g <gc_file> -p output_prefix [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=cooltoolsAB.__name__,
                        description=cooltoolsAB.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('-c','--cool',required=True,
                help='Input cool file')
    pReq.add_argument('-g','--gc',required=True,
                help='Input genome gc file')
    pReq.add_argument('-p','--prefix',required=True,
                      help='The prefix of the output file')
    pOpt.add_argument('-d','--output_Dir',default='.',
                      help='The output directory of the output file')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    check_file_exists(args.cool)
    check_file_exists(args.gc)
    log.info('------------ Start to call compartment ------------')
    clr = read_cool(args.cool)
    gc_cov = readGCfile(args.gc)
    log.info('Read GC content file done!')
    view_df = filter_chrom(clr)
    log.info('Filter chromosome done!')
    eigenvector_track = callCisEigs(clr,gc_cov,view_df)
    log.info('Call eigenvector done!')
    cvd = calExpectedCis(clr,view_df)
    log.info('Calculate expected done!')
    interaction_sum, interaction_count = calSaddle(clr,cvd,eigenvector_track,view_df)
    log.info('Calculate saddle done!')
    outputEigenvector(eigenvector_track,args.prefix,args.output_Dir)
    outputSum_Count(interaction_sum,interaction_count,args.prefix,args.output_Dir)


def plotSaddle(args):
    """
    Plot saddle plot based on cooltools results
    >>> %(prog)s -s <sum_file> -c <count_file> -e <eigenvector_file> -p output_prefix [Options]
    """
    install()
    p = argparse.ArgumentParser(prog=plotSaddle.__name__,
                        description=plotSaddle.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('-s','--sum',required=True,
                help='Input sum file')
    pReq.add_argument('-c','--count',required=True,
                help='Input count file')
    pReq.add_argument('-e','--eigenvector',required=True,
                help='Input eigenvector track file')
    pReq.add_argument('-p','--prefix',required=True,
                      help='The prefix of the output file')
    pOpt.add_argument('-d','--output_Dir',default='.',
                      help='The output directory of the output file')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    check_file_exists(args.sum)
    check_file_exists(args.count)
    check_file_exists(args.eigenvector)
    interaction_sum = readInteractionSum(args.sum)
    interaction_count = readInteractionCount(args.count)
    eigenvector_track = readEigenvector(args.eigenvector)
    plotsaddle(interaction_sum,interaction_count,eigenvector_track,args.prefix,args.output_Dir)


def ABStrength(args):
    """
    Calculate the Compartmentalisation strength based on cooltools results
    >>> %(prog)s -s <sum_file> -c <count_file> -p output_prefix [Options]
    """
    install()
    p = argparse.ArgumentParser(prog=ABStrength.__name__,
                        description=ABStrength.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('-s','--sum',required=True,
                help='Input sum file')
    pReq.add_argument('-c','--count',required=True,
                help='Input count file')
    pReq.add_argument('-p','--prefix',required=True,
                      help='The prefix of the output file')
    pOpt.add_argument('-d','--output_Dir',default='.',
                      help='The output directory of the output file')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    check_file_exists(args.sum)
    check_file_exists(args.count)
    interaction_sum = readInteractionSum(args.sum)
    interaction_count = readInteractionCount(args.count)
    strength = calSaddleStrength(interaction_sum,interaction_count)
    outputSaddleStrength(strength,args.prefix,args.output_Dir)


def main():
    actions = (
            ("CountsGC", "Counting GC content of each bin based on the bin of Cool and fasta file"),
            ("cooltoolsAB", "Call compartment based on cooltools algorithm"),
            ("plotSaddle", "Plot saddle plot based on cooltools results"),
            ("ABStrength", "Calculate the Compartmentalisation strength based on cooltools results"),
            ("CalBinFPKM", "Calculate the average expression value of each bin"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()