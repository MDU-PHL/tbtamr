## pip (simple install)

If you already have annotated vcf files and you do not need to do any variant calling or annotation, you can simply install the basic `tbtAMR` python package. **Just to be clear this is a bare bones installation and you will not be able to use fastq files or unannotated vcf files as inputs.**. If you require annotation - please see below.

You will need have python=3.10 or greater installed on your system (although management with mamba or conda is also ofcourse possible - actually recommended)

Optional:
```
mamba create -n tbtamr python=3.10
mamba activate tbtamr
```
OR 
```
conda create -n tbtamr python=3.10
conda activate tbtamr
```
Required:
```
pip3 install tbtamr
```


## Bioconda (recommended)

You can use `conda` to install the full `tbtAMR` package. This installation allows you to generate a vcf file from paired end reads, annotate a vcf file and/or generate genomic DST results.


```
conda (or mamba) create -n tbtamr tbtamr
conda activate tbtamr 
tbtamr -v
> 1.0.3
```


## Custom

### `tbtAMR` + annotation

In some cases you may have a vcf file - but need annotation (and do not want any of the variant calling dependencies). 

```
conda(mamba) create -n tbtamr snpEff
conda(mamba) activate tbtamr
pip3 install tbtamr==1.0.3
```


### `tbtamr` + lineage calling

If you would like to run lineage calling you can install [`pathogenprofiler`](https://github.com/jodyphelan/pathogen-profiler/tree/v4.3.0/pathogenprofiler) into your environment - note that this requires `bcftools`.

```
conda activate tbtamr
conda install bcftools (if you use conda to install tbtamr there is NO NEED for this step - only install bcftools if you use the pip3 simple install approach)
pip3 install git+https://github.com/jodyphelan/pathogen-profiler@v4.3.0
pip3 install pysam joblib tqdm pydantic requests
```

I have not suggested to install this with conda (which it absolutely can be if that works for you) as some of the dependency versions may clash with the `mutAMR` compatible dependencies. **If you wish to use `conda` to install `pathogen-profiler` please be aware that there may be unexpected behaviour**

