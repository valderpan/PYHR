# PYHR

**A custom genomics and epigenomics package**

## Function Overview
![PYHR_function](https://github.com/valderpan/PYHR/blob/master/PYHR.function.png)

## Dependencies

### requirement
- Python packages
    - natsort
    - intervaltree
    - argparse
    - path
    - pandas
    - numpy
    - biopython
    - gffutils
    - gzip
    - rich
    - scipy
- Executable programs
    - ice
    - build_matrix [from HiC-Pro/scripts/]

## Changelog
- v0.3.1: 
    - mv [3D.matrix] to [hic.matrix]
    - New module [3D.TAD]; new function [3D.TAD.boundaries2bed]
    - rename [hic.matrix CreateMatrix] to [hic.matrix CreateNNMatrix]
    - rename [hic.utils hicprostats] to [hic.utils statHiCPro]
