Common errors you may encounter
===============================

For errors that you don't find on this page or in the ``Snakemake`` (`documentation <https://snakemake.readthedocs.io/en/stable/>`_), we encourage you to submit a ticket through `GitHub issues <https://github.com/AlexanderLabWHOI/eukrhythmic/issues>`_. Please read the past issues first to see if we've addressed it before. If an open issue already exists, but it hasn't been answered, please feel free to add additional details to the thread. We will get back to you as soon as we can!

If you receive this error::
    
    Error: Directory cannot be locked. Please make sure that no other Snakemake process is trying to create the same files in the following directory:
    
    path/to/eukrhythmic/dir/eukrhythmic
    
    If you are sure that no other instances of snakemake are running on this directory, the remaining lock was likely caused by a kill signal or a power loss. It can be removed with the --unlock argument.
    
This means that ``eukrhythmic`` was run at one point and was not able to exit gracefully. Just run ``snakemake -s eukrhythmic --unlock`` to remove the lock on your directory.