Installing eukrhythmic
======================

`eukrhythmic` is designed as a modular pipeline that you can fully customize as the user of the resource. As such, we recommend setting up the base environment provided in the `environment.yaml` file prior to running anything, but this isn't strictly necessary; what you need is `Snakemake`, preferably newer than version 6, `Python`, preferably newer than version 3.8, some version of `pandas` and some version of `pyyaml`. Having these things installed will help later steps go as smoothly as possible.

Setting up a `conda` environment for running `Snakemake`
--------------------------------------------------------

Initialize the pipeline by setting up a ``conda`` environment, such that all the requested packages are loaded.::

    conda env create eukrhythmic --file environment.yaml
    
Downloading eukrhythmic
-----------------------

You may download the pipeline by cloning it directly from ``GitHub``.::

    git clone https://github.com/AlexanderLabWHOI/eukrhythmic
    cd eukrhythmic
    
All ``eukrhythmic`` commands are run from the base directory.

If you want to receive further information or plan on executing from the command line, go ahead and run::

    alias eukrhythmic='./bin/eukrhythmic.sh'
    
Then, you can execute `eukrhythmic -h` to see that the software is present in the workspace, and refer to "Running the pipeline from the command line" for more information on using command-line arguments with `eukrhythmic`.