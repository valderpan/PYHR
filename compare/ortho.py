#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/05/10


import re
import sys
import argparse
import pandas as pd
from natsort import natsorted
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify, read_file
from PYHR.utils.fasta import Fa2dict


log = richlog()

def countfamilies(GeneCountsdf):
    col_names = [col for col in GeneCountsdf][1:-1]
    species2OG = {}
    for colnames in col_names:
        species_name = re.findall('([0-9a-zA-Z\_\.]+)\.pep',colnames)[0]
        log.info('Get species ID : {}'.format(species_name))
        families_OG = GeneCountsdf[GeneCountsdf[colnames] !=0].iloc[:,0].tolist()
        with open(species_name+'.group.txt','w') as f:
            for OG in families_OG:
                f.write(OG+'\n')
        species2OG[species_name] = families_OG
    return species2OG


def Out_dataframe(species2OG):
    speciesL = []
    for s in species2OG.keys():
        qdf = pd.DataFrame({s:species2OG[s]})
        speciesL.append(qdf)
    resdf = pd.concat(speciesL,axis=1)
    resdf.to_csv('Orthogroups.GeneCount2venn.tab',sep='\t',header=True,index=False)


def Venn2GeneID(venndf,Orthodf,speciesID,output=None):
    sp2OG = {}
    resID = []
    for col in [column for column in venndf]:
        sp = natsorted([col]) if not '|' in col else natsorted(col.split('|'))
        sp2OG[','.join(sp)] = [id for id in venndf.loc[:,col].tolist() if not isinstance(id,float)]
    
    qID = venndf[speciesID].unique().tolist()
    qdf = Orthodf[Orthodf['Orthogroup'].isin(qID)]
    rawID = qdf['{}.pep'.format(speciesID)].unique().tolist()
    
    for i in rawID:
        if "," in i:
            rawIDL = i.split(', ')
            for j in rawIDL:
                resID.append(j)
        else:
            resID.append(i)

    for i in resID:
        print(i,file=output)


def Venn2Genenum(venndf,GeneCdf):
    from rich import print
    OGD = {}
    col_name = [column for column in venndf]
    for col in col_name:
        qOG = [i for i in venndf[col].tolist() if pd.isnull(i) == False]
        qdf = GeneCdf[GeneCdf['Orthogroup'].isin(qOG)]
        if qdf.shape[0] == len(qOG):
            Gtotal = qdf['Total'].sum()
            OGD[col] = [len(qOG),Gtotal]
        else:
            log.error('[bold magenta]{} rows != len(OGids) !!! Please check it ![/bold magenta]'.format(col))
    for item in sorted(OGD.items(),key=lambda x:x[1][0],reverse=False):
        print('[bold green]{}[/bold green]\t[bold yellow]{}[/bold yellow]\t[bold red]{}[/bold red]'.format(item[0],OGD[item[0]][0],OGD[item[0]][1]))



#outside command 
def ortho2venn(args):
    '''
    Gene family clustering based on orthofinder results
    >>> %(prog)s <Orthogroups.GeneCount.tsv> [Orthogroups.GeneCount2venn.tab] [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=ortho2venn.__name__,
                        description=ortho2venn.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('genecount', 
            help="Input Genecount.tsv file")
    pOpt.add_argument('-h', '--help', action='help',
            help='Show help message and exit.')
    
    args = p.parse_args(args)
    check_file_exists(args.genecount)
    df = read_file(args.genecount)
    species2OG = countfamilies(df)
    Out_dataframe(species2OG)


def extractVenn2GeneID(args):
    '''
    Extraction of species-specific gene IDs in gene families
    >>> %(prog)s venn.tsv Orthogroups.tsv speciesID [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=extractVenn2GeneID.__name__,
                        description=extractVenn2GeneID.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('venn', 
            help="Input venn.tsv file")
    pReq.add_argument('ortho',
            help='Input Orthogroups.tsv file')
    pReq.add_argument('speciesID',
            help='Input the prefix of the pep file, e.g. LA Os Sspon')
    pOpt.add_argument('-o', '--out', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='Show help message and exit.')

    args = p.parse_args(args)

    check_file_exists(args.venn)
    check_file_exists(args.ortho)
    venndf = read_file(args.venn)
    orthodf = read_file(args.ortho)
    Venn2GeneID(venndf,orthodf,args.speciesID,args.out)


def countGenenum4Venn(args):
    '''
    Count the number of genes in each venn subset
    >>> %(prog)s venn.tsv Orthogroups.GeneCount.tsv [Options]
    '''
    install()
    p = argparse.ArgumentParser(prog=countGenenum4Venn.__name__,
                        description=countGenenum4Venn.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('venn', 
            help="Input venn.tsv file")
    pReq.add_argument('GeneCount', 
            help="Input Orthogroups.GeneCount.tsv file")
    pOpt.add_argument('-h', '--help', action='help',
            help='Show help message and exit.')

    args = p.parse_args(args)

    check_file_exists(args.venn)
    check_file_exists(args.GeneCount)
    venndf = read_file(args.venn)
    GCdf = read_file(args.GeneCount)
    Venn2Genenum(venndf,GCdf)


def main():
    actions = (
            ("ortho2venn", "Gene family clustering based on orthofinder results"),
            ("extractVenn2GeneID", "Extraction of species-specific gene IDs in gene families"),
            ("countGenenum4Venn", "Count the number of genes in each venn subset"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()