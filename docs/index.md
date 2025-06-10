# tbtAMR

![Python package](https://img.shields.io/badge/python-3.10-blue.svg)
![CircleCI](https://dl.circleci.com/status-badge/img/gh/MDU-PHL/tbtamr/tree/master.svg?style=svg)

**Welcome to the home of tbtAMR the newest member of the tAMR family of tools - cousin of abritAMR.**

## What tbtAMR is.

tbtAMR was developed to address the need for a data-driven and flexible solution to inferring DST for _M. tuberculosis_ in a clinical and public health laboratory. tbtAMR generates genomic DST results that are suitable for reporting in a clinical setting. In order to facillitate this users are able to supply at a minimum a vcf file (annotated with snpEff) and either leverage a default mutational catalougue and criteria. Or provide their own catalogues and intepretative criteria in csv format. This means that there is no need for a developer to update logic or databases, which makes updates and therefore reverification much more straightforward and user-firendly in a CPHL setting. In addition, with additional dependencies, if required DST can also be reported from paired-end reads.

![flowdiagram](img/tbtamr_flow.svg)

## What tbtAMR is NOT.

`tbtAMR` is not

* A database or catalogue of mutations (although it can take these as inputs)
* For any sequences that are NOT M. tuberculosis
* A tool to identify relationships between sequences or samples
* A report generator - tbtAMR provides results in a way that can be easily incorporated into LIMS.
For detailed instructions for installation, basic usage and use of your own custom criteria check out the links on the left

## A note about vcf formats and annotation

Ultimately the new tbtAMR was born out of a need to have lightweight tool that could very easily slip into existing pipelines with few dependencies and very little overhead in terms of maintenance. At its core - it does not matter what sequencing technology is used - as long as there is a vcf file that is annotated correctly (for detailed information about vcf specifications please look [here](http://samtools.github.io/hts-specs/) and [here](https://pcingola.github.io/SnpEff/snpeff/inputoutput/)).

By default the assumption is that the format of the vcf file is correct and that annotation has been undertaken by snpEff in order to be compatible with the format in which variants are described in WHO catalogue version 2. If you want to use a different annotation format you will need to make sure that the catalogue of mutations you use is formatted in the same way. There is more about this in the Customising your inputs section.
