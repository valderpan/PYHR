#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/04/30


"""
the base libraries of formats module
modified from https://github.com/wangyibin/TDGP/blob/master/formats/base.py
"""


import sys
# sys.path.append('D:\\vscode\\MyScripts\\FAFU-CGB')
from PYHR.apps.base import richlog


log = richlog()

class BaseFile():
    '''
    定义一个读取文件的log
    '''
    def __init__(self, filename):
        if filename:
            log.debug("Loading file `{}`".format(filename))


class LineFile(BaseFile, list):
    '''
    BaseFile类的子类，可用于生成读取多个文件log'''
    def __init__(self,filename, comment=None, load=False):
        super(LineFile, self).__init__(filename)
        if load:
            pass

