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
    - cooler
    - cooltools
    - bioframe
    - cytoolz
    - matplotlib
    - pysam
- Executable programs
    - ice
    - build_matrix [from HiC-Pro/scripts/]

## Changelog
- v0.3.4
    - New function [3D.AB.CountsGC]
    - New function [3D.AB.cooltoolsAB]
    - New function [3D.AB.plotSaddle]
    - New function [3D.AB.ABStremgth]
    - DEBUG [3D.TAD boundaries2bed]
    - DEBUG [utils.kit CheckMd5]

- v0.3.3
    - Update [hic.utils statHiCPro] and [utils.kit CheckMd5]
    - Debug [3D.TAD.boundaries2bed]

- v0.3.2
    - move [rna.utils StatHisat2MappedRate] to [utils.bam StatHisat2MappedRate]
    - rename [utils.bam StatHisat2MappedRate] to [utils.bam StatMappedRate]
    - DEBUG [utils.bam StatMappedRate]

- v0.3.1: 
    - move [3D.matrix] to [hic.matrix]
    - New module [3D.TAD]; new function [3D.TAD.boundaries2bed]
    - rename [hic.matrix CreateMatrix] to [hic.matrix CreateNNMatrix]
    - rename [hic.utils hicprostats] to [hic.utils statHiCPro]
