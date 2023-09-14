# tbtAMR

`tbtAMR` implements TB-profiler and custom logic developed at MDU to identify mutations linked to AMR mechanisms in _M. tuberculosis_ and generate reports suitable for public health in Victoria. It may also be suitable for use in research settings.

**`tbtAMR` is now accredited to ISO15189 standard by NATA for use in Victoria Australia.** 


## tbtAMR installation

In order to install `tbtAMR` conda is strongly recommended - installation instructions can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

```
git clone git@github.com:kristyhoran/tbtamr.git
cd tbtamr
conda env create -f environment.yml
conda activate tbtamr
pip3 install .
tbtamr setup
tbtamr check
```

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
conda activate tb-profile-env
tb-profiler create_db -p <db_prefix> -c <path_to_csv> --custom --db_name <name_of_db> --db_commit <some_unique_string> --db_author <name_of_author> --db_date <date_created>
```

