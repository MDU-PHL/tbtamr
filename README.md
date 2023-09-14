[![CircleCI](https://dl.circleci.com/status-badge/img/gh/MDU-PHL/tbtamr/tree/master.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/MDU-PHL/tbtamr/tree/master)
[![Python 3.9](https://img.shields.io/badge/python-3.8-blue.svg)](https://www.python.org/downloads/release/python-390/)


# tbtAMR  <img src="https://github.com/MDU-PHL/tbtamr/blob/master/tbtamr_logo_transparent.png" width="100" height="70">

`tbtAMR` implements TB-profiler and custom logic developed at MDU to identify mutations linked to AMR mechanisms in _M. tuberculosis_ and generate reports suitable for public health in Victoria. It may also be suitable for use in research settings.

**`tbtAMR` is now accredited to ISO15189 standard by NATA for use in Victoria Australia.** 


## tbtAMR installation

In order to install `tbtAMR` conda is strongly recommended - installation instructions can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

```
conda create -n tbtamr tbtamr
conda activate tbtamr
tbtamr setup
tbtamr check
```

## Usage

At the moment `tbtAMR` is specific for illumina paired end reads. If you wish to use ONT data please use [TB-profiler](https://github.com/jodyphelan/TBProfiler). 

`tbtAMR` can be run with a single sample or in a batch approach. 

### Single sample

```
tbtamr run -r1 /path/to/R1.fq.gz -r2 /path/to/R2.fq.gz -px sample_name
```

### Batch 

```
tbtamr run -i /path/to/input.txt
```

This input file should have NO HEADER

## tbtAMR DB

`tbtAMR` comes with a modified mutational database, defined by validation at MDU for the purposes of reporting mutations in a public health and clinical setting in Victoria.

This database was created using `tb-profiler create_db` (from MDU fork of pathogen-profiler)

The following files are required in order to generate a custom database

* `genome.fasta` the reference for mutations in the database
* `genome.fasta.fai`
* `genome.gff`
* `barcode.bed` 
* a `csv` file
    * header ```Gene,Mutation,Drug,Confers,Interaction```

Creation can be done by navigating to a storage directory and running

```
conda activate tbtamr
tb-profiler create_db -p <db_prefix> -c <path_to_csv> --custom --db_name <name_of_db> --db_commit <some_unique_string> --db_author <name_of_author> --db_date <date_created>
```

