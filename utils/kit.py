#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/09/08


import sys
import argparse
from rich.traceback import install
from PYHR.apps.base import ActionDispatcher
from PYHR.apps.base import check_file_exists, richlog
from PYHR.apps.base import listify, read_file

log = richlog()

class md5():
    def __init__(self,file) -> None:
        self.md5D = {}
        self.file = file
    
    def set_md5D(self):
        with open(self.file) as f:
            lines = (line.strip() for line in f)
            for line in lines:
                md5_value,file_name = line.split()
                if not 'Rawdata' in file_name:
                    if '/' in file_name:
                        file_name = file_name.split('/')[-1]
                    self.md5D[file_name] = md5_value
        return self.md5D


def compare_md5(Omd5D,Nmd5D):
    Error_match = 0
    Correct_match = 0
    if len(Omd5D) != len(Nmd5D):
        log.error('The number of files before and after the transfer is different, please check it ！')
    else:
        log.debug('Consistent number of files before and after transfer ～')
        for key in Nmd5D.keys():
            if not key in Omd5D.keys():
                log.error('{} md5 value is missing'.format(key))
            if Nmd5D[key] == Omd5D[key]:
                log.info('{} md5 value is checked correctly'.format(key))
                Correct_match += 1
            else:
                print('='*30)
                log.error('{} md5 value is incorrectly checked'.format(key))
                log.debug('The original md5 value of file {} is {}'.format(key,Omd5D[key]))
                log.debug('The transferred md5 value of file {} is {}'.format(key,Nmd5D[key]))
                print('='*30)
                Error_match += 1
    log.info('Total file number : {},Correct checked file number : {}, Incorrect checked file number is {}'.format(len(Omd5D.keys()),
            Correct_match,Error_match))


## outside command 
def CheckMd5(args):
    """
    Correct the md5 value before and after file transfer
    >>> %(prog)s <Original md5 file> <transferred md5 file> [Options]
    """ 
    install()
    p = argparse.ArgumentParser(prog=CheckMd5.__name__,
                        description=CheckMd5.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('originalfile',
            help='Input original md5 file')
    pReq.add_argument('transferredfile',
            help='Input the transferred file')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    check_file_exists(args.originalfile)
    check_file_exists(args.transferredfile)
    od = md5(args.originalfile)
    nd = md5(args.transferredfile)
    compare_md5(od.set_md5D(),od.set_md5D())


def main():
    actions = (
            ("CheckMd5", "Correct the md5 value before and after file transfer"),
            #("concatKs2ggplotdensity", "Extract query seq by header"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()