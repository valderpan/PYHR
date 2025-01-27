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
    - coolpuppy
- Executable programs
    - ice
    - build_matrix [from HiC-Pro/scripts/]

## Changelog
- v0.5.1
    - DEBUG [utils.kit compare_md5]
    
- v0.4.6
    - DEBUG [utils.kit compare_md5] 
    - Rename [hic.utils] to [hic.base]
    - Rename [utils.Bam StatMappedRate] to [utils.Bam StatMappedRateHisat2]

- v0.4.5
    - DEBUG [threeD.loop read_loop AND findOverlapLoops]
    - Rename [utils.kit CheckMd5] to [utils.kit CheckMD5]

- v0.4.4
    - New function [hic.matrix Dump2Matrix]

- v0.4.3
    - DEBUG [utils.kit CheckMd5]
    - DEBUG [hic.loop findOverlapLoops]
    - New function [hic.utils DownsampleValid]

- v0.4.1
    - Rename [PYHR.3D] to [PYHR.threeD]
    - New module [threeD.base]
    - New function [threeD.TAD ATA]
    - Move [threeD.AB read_cool] to [threeD.base read_cool]
    - Move [threeD.AB readGCfile] to [threeD.base readGCfile]
    - Move [threeD.AB getResolution] to [threeD.base getResolution]
    - Move [threeD.AB countGC] to [threeD.base countGC]
    - Move [threeD.AB filter_chrom] to [threeD.base filter_chrom]
    - Move [threeD.AB calExpectedCis] to [threeD.base calExpectedCis]

- v0.3.4
    - New function [3D.AB CountsGC]
    - New function [3D.AB cooltoolsAB]
    - New function [3D.AB plotSaddle]
    - New function [3D.AB ABStremgth]
    - DEBUG [3D.TAD boundaries2bed]
    - DEBUG [utils.kit CheckMd5]

- v0.3.3
    - Update [hic.utils statHiCPro] and [utils.kit CheckMd5]
    - Debug [3D.TAD boundaries2bed]

- v0.3.2
    - move [rna.utils StatHisat2MappedRate] to [utils.bam StatHisat2MappedRate]
    - rename [utils.bam StatHisat2MappedRate] to [utils.bam StatMappedRate]
    - DEBUG [utils.bam StatMappedRate]

- v0.3.1: 
    - move [3D.matrix] to [hic.matrix]
    - New module [3D.TAD]; new function [3D.TAD boundaries2bed]
    - rename [hic.matrix CreateMatrix] to [hic.matrix CreateNNMatrix]
    - rename [hic.utils hicprostats] to [hic.utils statHiCPro]
