[![CircleCI](https://dl.circleci.com/status-badge/img/gh/MDU-PHL/tbtamr/tree/master.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/MDU-PHL/tbtamr/tree/master)
[![Python 3.x](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![ISO Accredited](https://img.shields.io/badge/ISO-15189_Accredited-green)](https://www.nata.com.au/)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/tbtamr.svg?label=Bioconda)](https://bioconda.github.io/recipes/tbtamr/README.html)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.landig.2025.100939-blue)](https://doi.org/10.1016/j.landig.2025.100939)

# tbtAMR  

<img src="https://github.com/MDU-PHL/tbtamr/blob/master/tbtamr_logo_transparent.png" align="right" width="100" height="70">

`tbtAMR` identifies and reports
mutations linked to anti-microbial resistance
in _M.tuberculosis_ from 
Illumina whole genome sequencing data.
It has been accredited to 
[ISO-15189 standards](https://en.wikipedia.org/wiki/ISO_15189)
by 
[NATA](https://nata.com.au/)
at the 
[MDU-PHL](https://biomedicalsciences.unimelb.edu.au/departments/microbiology-Immunology/research/services/microbiological-diagnostic-unit-public-health-laboratory).

## Installation

```
% conda create -n tbtamr -c bioconda tbtamr
% conda activate tbtamr
% tbtamr --version
```

## Quick start

```
% tbtamr full -s NAME -1 R1.fq.gz -2 R2.fq.gz
% cat btamr_results.csv

```

## Documentation

See our [wiki](https://github.com/MDU-PHL/tbtamr/wiki) page for further information on `tbtamr` usage.

## Feedback

File questions, bugs, or ideas on the 
[Issues page](https://github.com/MDU-PHL/tbtamr/issues).

## Licence

[GPLv3](https://raw.githubusercontent.com/MDU-PHL/tbtmar/master/LICENSE)

## Citation

Horan KA, Viberg L, Ballard SA, Globan M, 
Wirth W, Bond K, Webb JR, Dorji T, 
Williamson DA, Sait ML, Tay EL, 
Denholm JT, Howden BP, Seemann T, Sherry NL. 
_Bringing tuberculosis genomics to the clinic: development and validation of a comprehensive pipeline to predict antimicrobial susceptibility from genomic data, accredited to ISO standards._
**Lancet Digital Health.**
2025 Dec;7(12):100939.
PMID:[41436327](https://pubmed.ncbi.nlm.nih.gov/41436327/)

## Maintainer

[Krsity Horan](https://github.com/kristyhoran)
