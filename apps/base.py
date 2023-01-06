#!/usr/bin/env python3
# -*- encoding:utf-8 -*- 
# @Author : Haoran Pan
# Date: 2022/03/07


"""
Basic support for running library as script
"""


import os
import sys
import logging
import subprocess
import pandas as pd
import os.path as op
from natsort import natsorted
from PYHR import __copyright__, __version__
from rich.logging import Console, RichHandler


PYHRHELP = "PYHR utility libraries {} [{}]\n".format(__version__, __copyright__)


class ActionDispatcher(object):
    """
    This class will be invoked
    a) when the base package is run via __main__, listing all MODULESs
    a) when a directory is run via __main__, listing all SCRIPTs
    b) when a script is run directly, listing all ACTIONs
    This is controlled through the meta variable, which is automatically
    determined in get_meta().

    Copy from jcvi:https://github.com/tanghaibao/jcvi/blob/main/jcvi/apps/base.py
    """

    def __init__(self, actions):

        self.actions = actions
        if not actions:
            actions = [(None, None)]
        self.valid_actions, self.action_helps = zip(*actions)

    def get_meta(self):
        args = splitall(sys.argv[0])[-3:]
        args[-1] = args[-1].replace(".py", "")
        if args[-2] == "PYHR":
            meta = "MODULE"
        elif args[-1] == "__main__":
            meta = "SCRIPT"
        else:
            meta = "ACTION"
        return meta, args

    def print_help(self):
        meta, args = self.get_meta()
        if meta == "MODULE":
            del args[0]
            args[-1] = meta
        elif meta == "SCRIPT":
            args[-1] = meta
        else:
            args[-1] += " " + meta

        help = "Usage:\n    python -m {0}\n\n\n".format(".".join(args))
        help += "Available {0}s:\n".format(meta)
        max_action_len = max(len(action) for action, ah in self.actions)
        for action, action_help in sorted(self.actions):
            action = action.rjust(max_action_len + 4)
            help += (
                " | ".join((action, action_help[0].upper() + action_help[1:])) + "\n"
            )
        help += "\n" + PYHRHELP

        sys.stderr.write(help)
        sys.exit(1)

    def dispatch(self, globals):
        from difflib import get_close_matches

        meta = "ACTION"  # function is only invoked for listing ACTIONs
        if len(sys.argv) == 1:
            self.print_help()

        action = sys.argv[1]

        if not action in self.valid_actions:
            print("[error] {0} not a valid {1}\n".format(action, meta), file=sys.stderr)
            alt = get_close_matches(action, self.valid_actions)
            print(
                "Did you mean one of these?\n\t{0}\n".format(", ".join(alt)),
                file=sys.stderr,
            )
            self.print_help()

        globals[action](sys.argv[2:])


def richlog(level=logging.NOTSET):
    """
    Turn on the debugging
    """
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%Y-%m-%d %H:%M:%S]",
        handlers=[RichHandler(console=Console(stderr=True),rich_tracebacks=True)]
    )
    log = logging.getLogger("rich")
    return log

log = richlog()


def main():
    actions = (
        ('richlog', 'richlogsss'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def splitall(path):
    allparts = []
    while True:
        path, p1 = op.split(path)
        if not p1:
            break
        allparts.append(p1)
    allparts = allparts[::-1]
    return allparts


# def main():

#     actions = (
#             ("test", "test"),
#             ("plotLineRegress", "Plot two species pca1 lineregress"),
#             ('plotMultiLineRegress', "plot two species pca1 lineregress per switch type"),
#             ('plotLRPerChrom', 'plot two samples pca1 linregress per chromosome'),
#             ('plotEnrichment', 'plot compartment strength enrichment'),
#             ('plotStrength', 'plot compartment strength in multi samples'),
#             ('plotSwitchPie', 'plot two samples switch type pie picture'),
#             ('plotBoxPerChrom', 'plot boxplot of some data per chromosomes'),
#             ('plotBoxMultiSamples', 'plot boxplot of some data per samples on one chromosome'),
#             ('quickPlot', 'quick plot A/B compartments to show pc1 and gene density'),
#             ('getSyntenyGenePca', "get the synteny gene pairs pca value"),
#             ('annotateSwitchType', "annotate swithch type for synteny gene pairs"),
#             ('annotateType', 'annotate the compartment type for a pca bg file'),
#             ('statAB', 'stat A/B compartments informations'),
#             ('buscoGeneDist', 'to analysis busco gene distribution between A and B'),
#             ('getSwitchLink', 'to get links between sample1 and sample2 to plot circos')
            
#         )
#     p = ActionDispatcher(actions)
#     p.dispatch(globals())


def natglob(pathname, pattern=None):
    """
    Wraps around glob.glob(), but return a sorted list.
    """
    import glob as gl

    if pattern:
        pathname = op.join(pathname, pattern)
    return natsorted(gl.glob(pathname))


def get_module_docstring(filepath): #TODO
    """Get module-level docstring of Python module at filepath, e.g. 'path/to/file.py'."""
    co = compile(open(filepath).read(), filepath, "exec")
    if co.co_consts and isinstance(co.co_consts[0], str):
        docstring = co.co_consts[0]
    else:
        docstring = None
    return docstring


def dmain(mainfile, type="action"): #TODO
    cwd = op.dirname(mainfile)
    pyscripts = (
        [x for x in natglob(op.join(cwd, "*", "__main__.py"))]
        if type == "module"
        else natglob(op.join(cwd, "*.py"))
    )
    actions = []
    for ps in sorted(pyscripts):
        action = (
            op.basename(op.dirname(ps))
            if type == "module"
            else op.basename(ps).replace(".py", "")
        )
        if action[0] == "_":  # hidden namespace
            continue
        pd = get_module_docstring(ps)
        action_help = (
            [
                x.rstrip(":.,\n")
                for x in pd.splitlines(True)
                if len(x.strip()) > 10 and x[0] != "%"
            ][0]
            if pd
            else "no docstring found"
        )
        actions.append((action, action_help))

    a = ActionDispatcher(actions)
    a.print_help()



def listify(item):
    """
    To return a list or tuple value.
    将不是列表和元组的数据转换为列表
    From https://github.com/tanghaibao/jcvi/blob/master/jcvi/apps/base.py listify
    """
    return item if (isinstance(item, list) or isinstance(item, tuple)) else [item]


def check_file_exists(infile, ifcontinue=False):
    """
    Check whether the file or directory exists
    """
    if not op.exists(infile):
        if not ifcontinue:
            log.error('No such file of `{}`'.format(infile))
            sys.exit()
        else:
            log.warning('No such file of `{}`, but will continue'.format(infile))
        return False
    else:
        log.debug('Read file of `{}`'.format(infile))
        return True


def check_path_exists(path,ifcontinue=False):
    """
    check whether the path exists
    """
    if not os.path.exists(path):
        if ifcontinue == False:
            log.error('Path does not exist. Input path : {}'.format(path))
            sys.exit
        else:
            log.warning('Path does not exist, but will continue run. Input path : {}'.format(path))
        return False
    else:
        log.debug('Path exists, read path {}'.format(path))
        return True
        

def read_file(file,Names=None):
    if file.endswith('.xlsx') or file.endswith('.xls'):
        if Names:
            df = pd.read_excel(file,names=Names)
        else:
            df = pd.read_excel(file)
    elif file.endswith('.csv'):
        if Names:
            df = pd.read_csv(file,names=Names)
        else:
            df = pd.read_csv(file)
    else:
        if Names:
            df = pd.read_table(file,sep='\t',names=Names)
        else:
            df = pd.read_table(file,sep='\t')
    return df


def runshell(cmd):
    log.info(f'Run the shell command: {cmd}')
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        raise Exception(f"Shell command `{cmd}` failed to execute !")

if __name__ == "__main__":
    main()
    dmain()