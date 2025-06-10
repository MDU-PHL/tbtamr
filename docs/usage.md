## Before you start

`tbtAMR` can take any correctly formatted VCF file - which means that it is sequence platform agnostic. Alternatively, you can also supply Illumina paired-end reads as sqeuence input, providing you have installed `mutAMR` correctly. However when designing and developing `tbtAMR` some assumptions were made. 

- If you are using vcf as input - it is formatted correctly.
- Annotation was undertaken using snpEff or another annotation tool which adheres to specification found [here](https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf). If you supply an un-annotated vcf and have snpEff installed, `tbtAMR` will annotate for you (see HERE for installation instructions)

### Is your vcf file correctly formatted?
There is a standard for generation of VCF files, you can learn more about this [here](http://samtools.github.io/hts-specs/). Most bioinformatics tools will (should) follow this standard, if the tool you are using does not - there is no guarantee that `tbtAMR` behaviour will be as expected. Basic checks will be undertaken to ensure that your VCF is correct and if obvious issues are detected `tbtAMR` will error and ask you to check your inputs. However, if you are using a tool that does not follow this standard, but generates a file that mostly looks correct, `tbtAMR` may not detect the difference, it will appear to run correctly but the outputs may be unexpected or incorrect. If this is the case, please check your VCF generation and make sure that it follows the above linked specifications.

### Have you installed `tbtAMR` to run the way you want to run?

If you are supplying an annotated VCF as input and have no need for VCF generation or annotation, you can just install `tbtAMR` as describe [here](https://github.com/MDU-PHL/tbtamr/wiki/Installation#pip-simple-install). However, if you would like to supply paired-end fastq files as input, you will need to make sure that you have installed `tbtAMR` as descibed [here](https://github.com/MDU-PHL/tbtamr/wiki/Installation#bioconda-recommended).

Alternatively, you can also just install `snpEff` and `tbtAMR` and described [here](https://github.com/MDU-PHL/tbtamr/wiki/Installation#tbtamr--annotation) if you want `tbtAMR` to annotate your VCF file.

#### A note on lineage calling.

By default, `tbtAMR` does not call lineage for _M. tuberculosis_, you can do this if you like by following the installation instructions [here](https://github.com/MDU-PHL/tbtamr/wiki/Installation#tbtamr-and-lineage-calling) and using the `--call_lineage` flag.


## Generating genomic DST using default criteria (tbtAMR defaults)

### VCF input

```
tbtamr predict --vcf sample.vcf.gz -s sample_name
```

### paired-end reads input

```
tbtamr full -1 R1.fq.gz -2 R1.fq.gz -s sample_name
```
### Add phylogenetic lineage to outputs

Please make sure that you have installed `pathogen-profiler` as described [here](https://github.com/MDU-PHL/tbtamr/wiki/Installation#tbtamr-and-lineage-calling).

```
tbtamr predict --vcf sample.vcf.gz -s sample_name --call_lineage
```
OR

```
tbtamr full -1 R1.fq.gz -2 R1.fq.gz -s sample_name --call_lineage
```

## Generating genomic DST using custom criteria

Check out the instructions for defining criteria [here](https://github.com/MDU-PHL/tbtamr/wiki/Customising-your-inputs), example files can also be found [here](https://github.com/MDU-PHL/tbtamr/blob/master/example_criteria).

```
tbtamr predict --vcf sample.vcf.gz -s sample_name -c custom_mutational_catalogue.csv -cfg custom_db_config.json -r custom_interpretation_criteria.csv -cr custom_classification_criteria.csv
```
OR

```
tbtamr full -1 R1.fq.gz -2 R1.fq.gz -s sample_name -c custom_mutational_catalogue.csv -cfg custom_db_config.json -r custom_interpretation_criteria.csv -cr custom_classification_criteria.csv
```
## Use `tbtAMR` to generate annotated VCF file (no DST)

It is unclear why you may want to do this - but the functionality is here for you if you want to use it

```
tbtamr fq2vcf -1 R1.fq.gz -2 R1.fq.gz -s sample_name
```
OR
```
tbtamr annotate --vcf sample.vcf.gz -s sample_name
```

## Search your catalogue for variant information

`tbtAMR` also allows you to quickly peruse your catalogue - which can be useful if you don't want to break excel and you have a large catalogue.

```
tbtamr search -q rpoB_p.Ser540Leu rpoB_p.Ser540Ala
```
If you are using a custom catalogue - you will need to do 

```
tbtamr search -q some_query -c custom_mutational_catalogue.csv -cfg custom_db_config.json
```

## Output files

`tbtAMR` will output two files as standard. 

1. The first is `tbtamr_results.csv`. This is essentially the raw output of your interpretative criteria, this may include reportable and non-reportable results (depending on how the interpretative criteria and catalogue configuration are set up).

2. The second file is `tbtamr_linelist_report.csv`. This file is a filtered version of what is in the `tbtamr_results.csv`, with only results deemed appropriate for reporting (as defined by the interpretative criteria and catalogue configuration).

If you have added the `--cascade` flag, and have `cascade_reporting` set up in the catalogue configuration, a 3rd file will be created
3. `tbtamr_linelist_cascade_report.csv`, this is a further filtered version of the linelist report, with specific drugs reported based on the resistance profile detected. For more information see [here](https://github.com/MDU-PHL/tbtamr/wiki/Cascade-reporting) and [here](https://github.com/MDU-PHL/tbtamr/wiki/Customising-your-inputs#cascade-reporting).

Note all output files are designed to be simple **csv** to allow for straightforward management and incorporation into LIMS. `tbtAMR` does not make any assumptions as to how the reports should appear or be communicated.