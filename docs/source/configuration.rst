Advanced: Writing a configuration file
======================================

To write a configuration file for ``eukrhythmic``, you need to edit the ``config.yaml`` file included in the ``eukrhythmic`` base directory. This YAML-formatted file can be modified by changing the entries to the right of each of the colons in each line of the file.

.. _config: 

Configuration file entries
--------------------------

Below is a listing of each supported entry in the configuration file (``config.yaml`` in the base directory) and how to specify each flag when using the pipeline.


.. list-table:: Title
   :widths: 25 75
   :header-rows: 1
   
   * - Flag in file
     - Meaning & how to specify
   * - ``metaT_sample``
     - The name of the sample file containing sample ids to be used as unique identifiers in the pipeline, descriptive sample names, and input FASTA file names.
   * - ``inputDIR``
     - The file directory where the input data is found. Currently, should be specified with "/" separators, but no trailing "/". Should begin with "/" only if you are going to the root of your file system (not a  relative path).
   * - ``checkqual``
     - Boolean flag for whether to run quality checking with ``salmon``, ``QUAST``, ``BUSCO``, etc. on assemblies. If 1, these quality checks are performed.
   * - ``spikefile``
     - A path to a FASTA file containing the sequence of any spiking that might affect reads. This will depend on experimental setup. If the file is not valid (e.g., if this flag is set to 0), nothing is done.
   * - ``runbbmap``
     - A boolean flag to specify whether to use a spike file to drop spiked reads, according to what was done in your experiment. If 1,  the spikefile is used; otherwise, this  filtering is either not performed or is not used downstream in the pipeline  (depending on whether a spike file exists).
   * - ``kmers`` 
     - A list of *k*-mer sizes to use, where applicable, in assembly. These should all be integer values (default: 20, 50, 110). The median *k*-mer value in this list will be used when just 1 *k*-mer value is required.
   * - ``assemblers``
     - The assemblers to be used to assemble the metatranscriptomes (which will later be  merged). All of the specified assemblers in this list should have matching Snakemake rules in the ``modules`` folder of the main pipeline directory (named identically), as well as "clean" rules (explained below).
   * - ``jobname``
     - A descriptive name to be used to name jobs on your high-performance computing system, such that you can track the progress of your workflow.
   * - ``adapter``
     - Path to a FASTA file containing the adapter used during sequencing. Defaults to a  static adapter file in the ``static`` directory.
   * - ``separategroups``
     - A boolean flag. If 1, specified assembly  groups in the ``metaT_sample`` file are used to co-assemble raw files. Otherwise, each raw file is assembled separately regardless of what is specified in the "AssemblyGroup" column of the input file.
   * - ``outputDIR``
     - The path to a directory where all program output will be stored.
   * - ``assembledDIR``
     - The directory to move assembled files to, relative to the output directory. Defaults to "assembled"; not necessary to specify.
   * - ``renamedDIR``
     - The directory to move "renamed" files to  (which are files with the name of the  assembler added to each FASTA header), relative to the output directory. Defaults to "assembled"; not necessary to specify.
   * - ``scratch``
     - The location to move unnecessary intermediate files to after computation.