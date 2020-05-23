## Metatranscriptomic Pipelining

Documentation for `eukrhythmic` a metatranscriptomics analysis pipeline that enables users to run 

- Trimming
- Quality statistics
- Assembly (using several available assemblers)
- Assembly analysis / quality comparisons

## Initializing the pipeline

Initialize the pipeline by setting up a `conda` environment, such that all the requested packages are loaded. 

```
conda env create metatrans --file environment.yaml
```

## Pipeline schematic

![Schematic](/.infrastructure/eukrhythmic_revised.png)

## How to use the pipeline

To use the pipeline, the most important thing to do is to populate `config.yaml` with the paths to your particular input and output directories. Personalizing this will allow the pipeline to pull the relevant files when computing the results of individual rules, so this step is crucial. 

All input `fastq` files must be contained in the same directory, the `inputDIR` location in the `config.yaml` file. Only these metatranscriptomic data will be included in the analysis. These files do _not_, however, need to be located inside the `eukrhythmic` directory (and it is recommended that they are located elsewhere). 

The next thing that needs to be done is to produce the sample file, containing all of the relevant information to run the 

## File naming

Your file names should not subset one another. So if one file is called "supercool_sample", another should not be called "supercool_sample_2". However, if the two were called "supercool_sample_1" and "supercool_sample_2", this would be fine, because neither name is *entirely* found within the other.

## The `metaT_sample` file 

In the `config.yaml` file, there is a listing for a file called `metaT_sample` in the configuration file. This is essentially the data input source as far as what sample names you are expecting to include in your analysis, as well as any other information about the samples that you would like to be used. This is essential if you would like to apply groupings and co-assemble several samples together, and in general it is essential for the pipeline to work as intended. 

Depending on the application, some columns of this file must be added or may not be necessary. For example, for repeated samples in the same location, the latitude and longitude may not be necessary, because geographic variation in the metatranscriptomic assembly will not be evaluated. Any data that are not included in the default steps in the provided script can be excluded. 

As detailed in `scripts/make_sampleinfo.ipynb`, the minimum required to run the general (default) pipeline without any comparative analysis between samples are the "SampleName", "SampleID" and the "FastqFileNames". SampleName and SampleID can be identical if no such distinction exists in your sample. However, strictly, "SampleName" is a descriptive name for your samples, whereas "SampleID" is how the samples will be named throughout the pipeline, thus is is preferable to minimize special characters and length in the "SampleID" column, whereas "SampleName" may be more verbose.

### Autogeneration of full `metaT_sample` file

Using the script `scripts/autogenerate_metaT_sample.py`, you can autogenerate a working sample file for whatever files are present in the directory specified in your configuration file as `inputDIR`. This file should be run from the base `eukrhythmic` directory. At minimum, it requires a name for the output `metaT_sample` file as a parameter, which will be saved in the `input` directory, and then should be specified in the configuration file as the `metaT_sample` file. (TODO: in the future, autopopulate the config file with this entry and allow users to run the pipeline without even running this separately, with a default name for the `metaT_sample` file, as long as the `inputDIR` is specified). 

So within the base `eukrhythmic` directory, the following command may be run: 

```
python scripts/autogenerate_metaT_sample.py testsampledata.txt
```

Optionally, additional parameters may be provided. The second optional parameter is a file extension, which defaults to "fastq.gz". The third and fourth optional parameters are labels for forward and reverse reads (defaults to "\_1" and "\_2", respectively), and the fifth optional parameter is an additional file suffix used to split the filename (e.g. 001; defaults to the file extension; specifically important for single-end reads). 

### Autogeneration of "FastqFileNames" column with "SampleID" column
In /scripts/, there is a Python script called `make_sample_file.py` that will generate the `fastq` file names column, given your data input folder and an existing `metaT_sample` file that contains sample IDs. This is a good option if you only wish to run a subset of the files in your input data folder, and want to make sure the `fastq` file names are properly formatted.

If you use this script, all of the files with the fastq extension listed in your `INPUTDIR` that have a match to entries in your SampleID column will be included. Optionally, the script accepts up to two input arguments that specify how forward/reverse reads are labeled. By default, "\_R", "\_1", and "\_2" are searched for.

### Notes about manually creating `metaT_sample`
If you specify "FastqFileNames" manually, **ensure that the files are named uniquely and that the entire unique choice of name is specified in this column**. Filenames that match to more than two `fastq` files in your input directory will raise an exception.

In the event that you want more control over how your samples are named, use `scripts/make_sampleinfo.ipynb`. The AssemblyGroup column may be omitted if you do not mind if your samples are assigned assembly groups according to numbers 1-_n_, where _n_ is your number of samples, but in this case `separategroups` must be set to 0 in your `config.yaml` file. (_In the future, this may be updated to be done automatically if the column is absent, as well_). 

Once you have generated the `metaT_sample` file containing the information about your samples, and have populated `config.yaml` with the relevant directories, including, importantly, `outputDIR`, which will be the location of your results, you are ready to run the pipeline.

## Configuration file entries

Below is a listing of each supported entry in the configuration file (`config.yaml` in the base directory) and how to specify each flag when using the pipeline.

| Flag in file 	| Meaning & how to specify 	|
|------------------	|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|
| `metaT_sample` 	| The name of the sample file containing sample ids to be used as unique identifiers in the pipeline, descriptive sample names, and input FASTA file names. 	|
| `inputDIR` 	| The file directory where the input data is found. Currently, should be specified with "/" separators, but no trailing "/". Should begin with "/" only if you are going to the root of your file system (not a  relative path). 	|
| `checkqual` 	| Boolean flag for whether to run quality checking with `salmon`, `QUAST`, `BUSCO`, etc. on assemblies. If 1, these quality checks are performed. 	|
| `spikefile` 	| A path to a FASTA file containing the sequence of any spiking that might affect reads. This will depend on experimental setup. If the file is not valid (e.g., if this flag is set to 0), nothing is done. 	|
| `dropspike` 	| A boolean flag to specify whether to use a spike file to drop spiked reads, according to what was done in your experiment. If 1,  the spikefile is used; otherwise, this  filtering is either not performed or is not used downstream in the pipeline  (depending on whether a spike file exists). 	|
| `kmers` 	| A list of *k*-mer sizes to use, where applicable, in assembly. These should all be integer values (default: 20, 50, 110). 	|
| `assemblers` 	| The assemblers to be used to assemble the metatranscriptomes (which will later be  merged). All of the specified assemblers in this list should have matching Snakemake rules in the `modules` folder of the main pipeline directory (named identically), as well as "clean" rules (explained below). 	|
| `jobname` 	| A descriptive name to be used to name jobs on your high-performance computing system, such that you can track the progress of your workflow. 	|
| `adapter` 	| Path to a FASTA file containing the adapter used during sequencing. Defaults to a  static adapter file in the `static` directory. 	|
| `separategroups` 	| A boolean flag. If 1, specified assembly  groups in the `metaT_sample` file are used to co-assemble raw files. Otherwise, each raw file is assembled separately regardless of what is specified in the "AssemblyGroup" column of the input file. 	|
| `outputDIR` 	| The path to a directory where all program output will be stored. 	|
| `assembledDIR` 	| The directory to move assembled files to, relative to the output directory. Defaults to "assembled"; not necessary to specify. 	|
| `renamedDIR` 	| The directory to move "renamed" files to  (which are files with the name of the  assembler added to each FASTA header), relative to the output directory. Defaults to "assembled"; not necessary to specify. 	|
| `scratch` 	| The location to move unnecessary intermediate files to after computation. 	|

## Running the pipeline

Once the pieces are in place, and you have either activated an environment using `environment.yaml` or otherwise installed `snakemake`, you can run the pipeline using:

```
sbatch submit/snake_submit.sh
```

If you are using the `SLURM` scheduler, or by simply executing the `submit/snake_submit.sh` file in the `eukrhythmic` directory if you are not using a scheduler, or are logged into a computer with sufficient computational resources (e.g., a `SLURM` job run in interactive mode). 

## Adding new assemblers

In order to add a new assembler to the pipeline, three things need to be included:

1. A rule for the assembler, including a command that can be run using software in the specified `conda` environment, called `<assembler>`.
2. A rule to "clean" the assembly output, named `<assembler>_clean`, including moving the completed assemblies to the shared assembly folder, specified as `assembledDIR`, which is a subdirectory of the `outputDIR`, also specified in the configuration file. Intermediate files should also be moved to a scratch directory or deleted, based on user preferences and space requirements. Any other files needed by other tools or desired by the user should be moved to a subdirectory of the output directory. If they are specified as output files, `snakemake` will generate them automatically. Otherwise, the user will need to manually create directories that do not already exist (specifying them as output files is more extensible). 
3. A list entry for the assembler in the configuration file that matches the name of `<assembler>` in each `snakemake` rule. 
