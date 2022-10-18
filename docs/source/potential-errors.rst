Common errors you may encounter
===============================

If you receive this error::
    
    Error: Directory cannot be locked. Please make sure that no other Snakemake process is trying to create the same files in the following directory:
    
    path/to/eukrhythmic/dir/eukrhythmic
    
    If you are sure that no other instances of snakemake are running on this directory, the remaining lock was likely caused by a kill signal or a power loss. It can be removed with the --unlock argument.
    
This means that ``eukrhythmic`` was run at one point and was not able to exit gracefully. Just run ``snakemake -s eukrhythmic --unlock`` to remove the lock on your directory.