Intermediate files and cleanup
==============================

.. _intermediate: 

A "scratch" directory, specified in the configuration file, is used to store intermediate files after execution of rules such as assembly, which produce many files which are not needed downstream in the pipeline. To override this behavior, specify the output directory and the scratch directory to be the same location.

After the pipeline has been run, simply enter::

    snakemake hardclean --cores 1

To safely remove the scratch directory, if you don't need the intermediate files generated in individual pipeline steps.