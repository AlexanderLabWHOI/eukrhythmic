import io
import os
from os import listdir
from os.path import isfile, join
import sys
import pandas as pd
import numpy as np
import pathlib
import yaml
from snakemake.exceptions import print_exception, WorkflowError  
sys.path.insert(1, 'scripts') 
from importworkspace import *

## THIS MAY ALL BE DELETED IN FUTURE ##
with open('config.yaml') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

# Read relevant values from the configuration file
KMERVALS = list(config['kmers'])
ASSEMBLERS = list(config['assemblers'])
OUTPUTDIR = config["outputDIR"]
SCRATCHDIR = config["scratch"]
INPUTDIR = config["inputDIR"]
DATAFILE = config["metaT_sample"]
if 'renamedDIR' in config:
    RENAMEDDIR = os.path.join(OUTPUTDIR, config['renamedDIR'])
else:
    RENAMEDDIR = os.path.join(OUTPUTDIR, "renamed")
if 'assembledDIR' in config:
    ASSEMBLEDDIR = os.path.join(OUTPUTDIR, config['assembledDIR'])
else:
    ASSEMBLEDDIR = os.path.join(OUTPUTDIR, "assembled")
directories = [ASSEMBLEDDIR,INPUTDIR,OUTPUTDIR,SCRATCHDIR,RENAMEDDIR]
## END SECTION TO BE DELETED ##

def checkrequirementsfct(): 
    print("\033[1;35m Checking that configuration file entries are valid...  \n")
    # Check to make sure a list of integer k-mer values is supplied.
    if isinstance(KMERVALS, list):
        for kmer in KMERVALS:
            try:
                if not isinstance(kmer + 1, int):
                    raise ValueError
            except ValueError:
                print("Please only specify integer k-mer values.")
                sys.exit(1)

    # Check to make sure a valid list of assemblers is supplied.
    if isinstance(ASSEMBLERS, list):
        for assembler_ind in range(0,len(ASSEMBLERS)):
            assembler = ASSEMBLERS[assembler_ind]
            try:
                assemblerrules = os.listdir("modules/")
                if assembler + "-snake" not in assemblerrules:
                    if assembler.lower() + "-snake" in assemblerrules:
                        ASSEMBLERS[assembler_ind] = assembler.lower()
                    else:
                        raise ValueError
            except ValueError:
                print("Please only specify assemblers for which rules are available.")
                sys.exit(1)
    else:
        try:
            raise ValueError
        except ValueError:
            print("Please specify your assemblers as a list. See examples for details.")
            sys.exit(1)

    # Check to make sure the user hasn't added trailing /.
    for dir_curr in directories:
        try:
            if ((dir_curr[len(dir_curr)-1] == '/') | (dir_curr[len(dir_curr)-1] == '\\')):
                raise ValueError
        except ValueError:
            print("Please do not add trailing slashes to input, output, scratch, assembled, or renamed directories.")
            sys.exit(1)
    
    # Check whether the required number of entries of each fastq file can be found in the input directory.
    inputfiles = "|".join(os.listdir(INPUTDIR))
    filenames = []
    singleorpaired = []
    for currfile_ind in range(0, len(fastqnames)):
        currfile = fastqnames[currfile_ind]
        occurrences = inputfiles.count(currfile)
        if occurrences > 2:
            try:
                raise ValueError
            except ValueError:
                print("There are too many occurrences of fastq file " + currfile + " in input directory.")
                sys.exit(1)
               
    
    print("\033[1;32m All initial checks complete.  \n")