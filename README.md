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





