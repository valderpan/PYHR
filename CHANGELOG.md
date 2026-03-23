## Changelog

All notable changes to this project will be documented in this file.
## [0.5.7] - 2026-03-12
### Added
- New function `utils.kit.DownloadGSESuppl` to download the supplementary files of GSE from GEO database.

### Changed
- Change the output format of the `rna.exp.stringtie2ExpMatrix` function from `.xlsx` to `.csv`


## [0.5.6] - 2025-12-25
### Added
- Implemented colorized terminal output for `utils.kit.CheckMD5` to improve status visibility.

### Changed
- Updated documentation for `hic.matrix.ICEMatrix`.


## [0.5.5] - 2025-08-21
### Added
- New function `apps.base.remove_suffix` for handling file extensions.
- New function `hic.base.estimateHiCresolution` to automatically calculate optimal resolution.

### Changed
- **Performance:** Switched to `pfastq-dump` instead of standard `fastq-dump` in `utils.kit.DownloadFastq` to accelerate SRA downloads.

### Fixed
- Fixed parsing logic in `threeD.loop.parse_inter` and `run_common` causing crashes on empty inputs.


- 0.5.5
    - New function [apps.base remove_suffix]
    - New function [hic.base estimateHiCresolution]
    - Use pfastq-dump instead of fastq-dump in [utils.kit DownloadFastq]
- 0.5.3
    - Update [apps.base runshell]
    - New module [utils.vcf]
    - New function [utils.vcf replaceVCFsampleID]
- 0.5.2
    - DEBUG [threeD.loop parse_inter AND run_common]

- 0.5.1
    - DEBUG [utils.kit compare_md5]
    
- 0.4.6
    - DEBUG [utils.kit compare_md5] 
    - Rename [hic.utils] to [hic.base]
    - Rename [utils.Bam StatMappedRate] to [utils.Bam StatMappedRateHisat2]

- 0.4.5
    - DEBUG [threeD.loop read_loop AND findOverlapLoops]
    - Rename [utils.kit CheckMd5] to [utils.kit CheckMD5]

- 0.4.4
    - New function [hic.matrix Dump2Matrix]

- 0.4.3
    - DEBUG [utils.kit CheckMd5]
    - DEBUG [hic.loop findOverlapLoops]
    - New function [hic.utils DownsampleValid]

- 0.4.1
    - Rename [PYHR.3D] to [PYHR.threeD]
    - New module [threeD.base]
    - New function [threeD.TAD ATA]
    - Move [threeD.AB read_cool] to [threeD.base read_cool]
    - Move [threeD.AB readGCfile] to [threeD.base readGCfile]
    - Move [threeD.AB getResolution] to [threeD.base getResolution]
    - Move [threeD.AB countGC] to [threeD.base countGC]
    - Move [threeD.AB filter_chrom] to [threeD.base filter_chrom]
    - Move [threeD.AB calExpectedCis] to [threeD.base calExpectedCis]

- 0.3.4
    - New function [3D.AB CountsGC]
    - New function [3D.AB cooltoolsAB]
    - New function [3D.AB plotSaddle]
    - New function [3D.AB ABStremgth]
    - DEBUG [3D.TAD boundaries2bed]
    - DEBUG [utils.kit CheckMd5]

- 0.3.3
    - Update [hic.utils statHiCPro] and [utils.kit CheckMd5]
    - Debug [3D.TAD boundaries2bed]

- 0.3.2
    - move [rna.utils StatHisat2MappedRate] to [utils.bam StatHisat2MappedRate]
    - rename [utils.bam StatHisat2MappedRate] to [utils.bam StatMappedRate]
    - DEBUG [utils.bam StatMappedRate]

- 0.3.1: 
    - move [3D.matrix] to [hic.matrix]
    - New module [3D.TAD]; new function [3D.TAD boundaries2bed]
    - rename [hic.matrix CreateMatrix] to [hic.matrix CreateNNMatrix]
    - rename [hic.utils hicprostats] to [hic.utils statHiCPro]