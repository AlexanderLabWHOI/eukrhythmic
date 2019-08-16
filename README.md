## Metatranscriptomic Pipelining

Playing around with creating a metatranscriptomic pipeline (for now)

## Initializing the pipeline

Initialize the pipeline by setting up a `conda` environment, such that all the requested packages are loaded. 

```
conda env create metatrans --file environment.yaml
```

## How to use the pipeline

To use the pipeline, the most important thing to do is to populate `config.yaml` with the paths to your particular input and output directories. Personalizing this will allow the pipeline to pull the relevant files when computing the results of individual rules, so this step is crucial. 

All input `fastq` files must be contained in the same directory, the `inputDIR` location in the `config.yaml` file. Only these metatranscriptomic data will be included in the analysis. 

## The `metaT_sample` file 

In the `config.yaml` file, there is a listing for a file called `metaT_sample` in the configuration file. This is essentially the data input source as far as what sample names you are expecting to include in your analysis, as well as any other information about the samples that you would like to be used in the analysis. 

Depending on the application, some columns of this file must be added or may not be necessary. For example, for repeated samples in the same location, the latitude and longitude may not be necessary, because geographic variation in the metatranscriptomic assembly will not be evaluated. 

The minimum required to run the general (default) pipeline without any comparative analysis between samples are the "SampleID" and the "FastqFileNames". In /scripts/, there is a Python script that will generate the `fastq` file names column, given your data input folder. If you use this script, all of the files with the fastq extension listed in your `INPUTDIR` that have a match to entries in your SampleID column will be included. This script also has a toggle for creating the metaT\_sample file essentially from scratch, in the case that you (1) want all `fastq` files to be included and (2) do not care about the choice of SampleID. 

## Generating sample file with the file names

Executing the `make_sample_file.py` file from the main directory of the repository (e.g. `python scripts/make_sample_file.py` will fit your input data file with the fastq file names from the directory you have specified as the input directory in `config.yaml`. 


