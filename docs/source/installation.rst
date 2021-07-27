Installing eukrhythmic
======================

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