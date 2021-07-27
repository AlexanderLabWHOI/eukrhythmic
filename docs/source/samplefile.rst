Sample file specification
=========================

.. _naming: 

File naming
-----------

Your file names should not subset one another. So if one file is called "supercool_sample", another should not be called "supercool_sample_2". However, if the two were called "supercool_sample_1" and "supercool_sample_2", this would be fine, because neither name is *entirely* found within the other.

.. _sample:

The ``metaT_sample`` file
-------------------------

In the `config.yaml` file, there is a listing for a file called `metaT_sample` in the configuration file. This is essentially the data input source as far as what sample names you are expecting to include in your analysis, as well as any other information about the samples that you would like to be used. This is essential if you would like to apply groupings and co-assemble several samples together, and in general it is essential for the pipeline to work as intended. 

Depending on the application, some columns of this file must be added or may not be necessary. For example, for repeated samples in the same location, the latitude and longitude may not be necessary, because geographic variation in the metatranscriptomic assembly will not be evaluated. Any data that are not included in the default steps in the provided script can be excluded. 

As detailed in `scripts/make_sampleinfo.ipynb`, the minimum required to run the general (default) pipeline without any comparative analysis between samples are the "SampleName", "SampleID" and the "FastqFileNames". SampleName and SampleID can be identical if no such distinction exists in your sample. However, strictly, "SampleName" is a descriptive name for your samples, whereas "SampleID" is how the samples will be named throughout the pipeline, thus is is preferable to minimize special characters and length in the "SampleID" column, whereas "SampleName" may be more verbose.

.. _fullauto:

Autogeneration of full ``metaT_sample`` file
--------------------------------------------

Using the script ``scripts/autogenerate_metaT_sample.py``, you can autogenerate a working sample file for whatever files are present in the directory specified in your configuration file as ``inputDIR``. This file should be run from the base `eukrhythmic` directory. At minimum, it requires a name for the output ``metaT_sample`` file as a parameter, which will be saved in the `input` directory, and then should be specified in the configuration file as the `metaT_sample` file. (Note: in the future, autopopulate the config file with this entry and allow users to run the pipeline without even running this separately, with a default name for the ``metaT_sample`` file, as long as the ``inputDIR`` is specified). 

So within the base ``eukrhythmic`` directory, the following command may be run::

    python scripts/autogenerate_metaT_sample.py testsampledata.txt

Optionally, additional parameters may be provided. The second optional parameter is a file extension, which defaults to "fastq.gz". The third and fourth optional parameters are labels for forward and reverse reads (defaults to "\_1" and "\_2", respectively), and the fifth optional parameter is an additional file suffix used to split the filename (e.g. 001; defaults to the file extension; specifically important for single-end reads). 

.. _fastqauto:

Autogeneration of "FastqFileNames" column with "SampleID" column
----------------------------------------------------------------

In /scripts/, there is a Python script called ``make_sample_file.py`` that will generate the ``fastq`` file names column, given your data input folder and an existing ``metaT_sample`` file that contains sample IDs. This is a good option if you only wish to run a subset of the files in your input data folder, and want to make sure the ``fastq`` file names are properly formatted.

If you use this script, all of the files with the fastq extension listed in your ``INPUTDIR`` that have a match to entries in your SampleID column will be included. Optionally, the script accepts up to two input arguments that specify how forward/reverse reads are labeled. By default, "\_R", "\_1", and "\_2" are searched for.

.. _manual:

Notes about manually creating ``metaT_sample``
----------------------------------------------

If you specify "FastqFileNames" manually, **ensure that the files are named uniquely and that the entire unique choice of name is specified in this column**. Filenames that match to more than two ``fastq`` files in your input directory will raise an exception.

In the event that you want more control over how your samples are named, use ``scripts/make_sampleinfo.ipynb``. The AssemblyGroup column may be omitted if you do not mind if your samples are assigned assembly groups according to numbers 1-*n*, where *n* is your number of samples, but in this case ``separategroups`` must be set to 0 in your ``config.yaml`` file. (*In the future, this may be updated to be done automatically if the column is absent, as well*). 

Once you have generated the ``metaT_sample`` file containing the information about your samples, and have populated ``config.yaml`` with the relevant directories, including, importantly, ``outputDIR``, which will be the location of your results, you are ready to run the pipeline.