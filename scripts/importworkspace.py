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

print("\033[1;35m Reading in variables from configuration file...  \n")
with open('config.yaml') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    
required_entries = ["metaT_sample","inputDIR","checkqual","spikefile","dropspike",\
                    "kmers","assemblers","jobname","adapter","separategroups","outputDIR","scratch",\
                   "transdecodercutoff","secondclustercutoff"]
for r in required_entries:
    try:
        if r not in config:
            offender = r
            raise ValueError
    except ValueError:
        print("Check that you have specified values in the configuration file. The missing entry that triggered this error was "  + str(offender) + ".")
        sys.exit(1)
        
## READ ALL RELEVANT DATA IN FROM CONFIGURATION FILE ##
if "metaT_sample" in config:
    DATAFILE = config["metaT_sample"]
else:
    DATAFILE = os.path.join("input", "sampledata.txt")

INPUTDIR = config["inputDIR"]
OUTPUTDIR = config["outputDIR"]
SCRATCHDIR = config["scratch"]
KMERVALS = list(config['kmers'])
SPIKEFILE = config['spikefile']
DROPSPIKE = config['dropspike']
TRANSDECODERORFSIZE = config["transdecodercutoff"]
MINCOVERAGECLUST2 = config["transdecodercutoff"]
ADAPTER = config["adapter"]
if os.path.isfile(SPIKEFILE):
    ISFILESPIKE = 1
else:
    ISFILESPIKE = 0
ASSEMBLERS = [curr.lower() for curr in list(config['assemblers'])]
SAMPLEINFO = pd.read_csv(DATAFILE, sep = "\t")
ASSEMBLYDICT = dict(zip(list(SAMPLEINFO.AssemblyGroup), list(SAMPLEINFO.SampleName)))
samplenames = list(SAMPLEINFO.SampleID);
fastqnames = list(SAMPLEINFO.FastqFile);

# Check to make sure the scratch directory exists; otherwise, create it.
os.system("mkdir -p " + SCRATCHDIR)

if 'renamedDIR' in config:
    RENAMEDDIR = os.path.join(OUTPUTDIR, config['renamedDIR'])
else:
    RENAMEDDIR = os.path.join(OUTPUTDIR, "renamed")
if 'assembledDIR' in config:
    ASSEMBLEDDIR = os.path.join(OUTPUTDIR, config['assembledDIR'])
else:
    ASSEMBLEDDIR = os.path.join(OUTPUTDIR, "assembled")
    
if DROPSPIKE == 0:
    LEFTFILE = lambda filename: expand(os.path.join(OUTPUTDIR, "firsttrim", "{samples}_1.trimmed.fastq.gz"), samples = get_samples(filename.assembly))
    RIGHTFILE = lambda filename: expand(os.path.join(OUTPUTDIR, "firsttrim", "{samples}_2.trimmed.fastq.gz"), samples = get_samples(filename.assembly))
else:
    LEFTFILE = lambda filename: expand(os.path.join(OUTPUTDIR, "bbmap", "{samples}_1.clean.fastq.gz"), samples = get_samples(filename.assembly))
    RIGHTFILE = lambda filename: expand(os.path.join(OUTPUTDIR, "bbmap", "{samples}_2.clean.fastq.gz"), samples = get_samples(filename.assembly))

directories = [ASSEMBLEDDIR,INPUTDIR,OUTPUTDIR,SCRATCHDIR,RENAMEDDIR]
        
print("\033[1;35m Checking directory formatting...  \n")
# Check to make sure the user hasn't added trailing /.
for dir_curr in directories:
    try:
        if ((dir_curr[len(dir_curr)-1] == '/') | (dir_curr[len(dir_curr)-1] == '\\')):
            raise ValueError
    except ValueError:
        print("Please do not add trailing slashes to input, output, scratch, assembled, or renamed directories.")
        sys.exit(1)

        
print("\033[1;35m Checking that appropriate input files exist...  \n")
inputfiles = "|".join(os.listdir(INPUTDIR))
filenames = []
singleorpaired = []
for currfile_ind in range(0, len(fastqnames)):
    currfile = fastqnames[currfile_ind]
    occurrences = inputfiles.count(currfile)
    if occurrences > 2:
        print("There are too many occurrences of fastq file " + currfile + " in input directory.")
    elif occurrences == 2:
        singleorpaired.extend([1,2])
        filenames.extend([samplenames[currfile_ind]] * 2)
    else:
        singleorpaired.append(1)
        filenames.append(samplenames[currfile_ind])

if config["separategroups"] == 1:
    assemblygroups = list(set(SAMPLEINFO.AssemblyGroup))
else:
    assemblygroups = [1] * len(INPUTFILES)

def combineassemblers(assembly):
    return(" ".join([os.path.join(ASSEMBLEDDIR, assembly + "_" + curr + ".fasta") for curr in ASSEMBLERS]))
