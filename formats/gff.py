#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/03/09


'''
Util tools for gff formats.
'''


import sys
import gffutils
import argparse
import pandas as pd
import os.path as op
from rich.traceback import install
from collections import OrderedDict
from PYHR.apps.base import ActionDispatcher,richlog


log = richlog()


def import_gff(gff_path, _type=None):
    compression = 'gzip' if gff_path.endswith('.gz') else 'infer'
    gff_df = pd.read_csv(gff_path, sep='\t', 
                         header=None, 
                         index_col=0, 
                         comment="#", 
                         compression = compression,
                         names=('seqid', 'source',
                                'type', 'start', 
                                'end', 'score', 
                                'strand','phase',
                                'attributes'))
    if _type: 
        if _type in set(gff_df.head(100)['type']):
            gff_df = gff_df[gff_df['type'] == _type]
        else:
            log.warning('Warning: Failed to filter data by {}, input type is not correct'.format(_type))
    return gff_df 


def get_gene_id(attributes):
    """
    取出ID对应的geneID
    """
    db = dict(map(lambda x: x.split('='), 
                  [i for i in attributes.split(';') if i]))
    return db['ID']


def createDB(gff, dbfn=None):
    """
    create a db file for gffutils
    Params:
    --------
    gff: `str` gff file
    dbfn: `str` db name [default: gff.db]
    Returns:
    --------
    db: `gffutils.FeatureDB` db for gffutils
    Examples:
    ---------
    >>> createDB('in.gff3')
    <FeatureDB ...>
    """
    ## create gff db
    if not dbfn: 
        db_name = gff + ".db"
    else:
        db_name = dbfn
    if not op.exists(db_name):
        log.debug('No such database file of `{}`, creating ...'.format(db_name))
        gffutils.create_db(gff, dbfn=db_name, keep_order=True)
    else:
        log.debug('Already exists DB file of `{}`, skip.'.format(db_name))
    db = gffutils.FeatureDB(db_name)
    
    return db


## outside command 
def RenameAttributesID(args):
    """
    Change the attributes within the gff3
    >>> %(prog)s <in.gff3> <rename.list> [Options]
    """
    install()
    p = argparse.ArgumentParser(prog=RenameAttributesID.__name__,
                        description=RenameAttributesID.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('gff', 
            help='gff file ')
    pReq.add_argument('renameList', 
            help='rename list, two columns <old_gene_name\tnew_gene_name>')
    pOpt.add_argument('-d', '--dbname', default=None,
            help='gff database name [default: gff3.db')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    output = args.output
    db = createDB(args.gff, args.dbname)
    rename_db = OrderedDict(i.strip().split() for i in open(args.renameList)
                        if i.strip())
    
    for ID in rename_db:
        print(str(db[ID]).replace(ID, rename_db[ID]), file=output)
        for feature in db.children(ID, order_by='start'):
            print(str(feature).replace(ID, rename_db[ID]), file=output)
        print("", file=output)
    log.debug("Successful. Output file is in `{}`".format(output.name))


def main():
    actions = (
            ("RenameAttributesID", "rename the ID in attributes"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()