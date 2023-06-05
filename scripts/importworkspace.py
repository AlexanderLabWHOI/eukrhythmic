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
from snakemake.logging import logger

#arguments = sys.argv
configfile = "config.yaml" #"/vortexfs1/omics/alexander/akrinos/2021-tara-phaeo/2021-akrinos-tara-phaeo/eukrhythmic_setup/config.yaml"
#if len(arguments) > 2:
#    configfile = arguments[2]
#    print(configfile,flush=True)

logger.info("\033[1;35m Reading in variables from configuration file...  \n")
with open(configfile) as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    
logger.info("\033[1;35m Reading in variables from cluster configuration file...  \n")
with open('cluster.yaml') as f:
    cluster = yaml.load(f, Loader=yaml.FullLoader)

required_entries = ["maxmemory","maxtasks","maxcores"]
if (cluster != None):
    if ("required" in cluster):
        for r in required_entries:
            try:
                if r not in cluster["required"].keys():
                    offender = r
                    raise ValueError
            except ValueError:
                logger.info("Check that you have specified values in the cluster configuration file file. The missing entry that triggered this error was "  + str(offender) + ".")
                sys.exit(1)

if "required" in cluster:
    DEFAULTQUEUE = cluster["required"]["defaultqueue"]
    DEFAULTUSER = cluster["required"]["accountname"]
    MAXMEMORY = cluster["required"]["maxmemory"]
    MAXTASKS = cluster["required"]["maxtasks"]
    MAXCORES = cluster["required"]["maxcores"]
    MAXCPUSPERTASK = cluster["required"]["maxcpuspertask"]
    
required_entries = ["metaT_sample","inputDIR","checkqual","spikefile","runbbmap",\
                    "kmers","assemblers","jobname","adapter","separategroups","outputDIR","scratch",\
                   "transdecodercutoff","secondclustercutoff","defaultkmer"]
for r in required_entries:
    try:
        if r not in config:
            offender = r
            raise ValueError
    except ValueError:
        logger.info("Check that you have specified values in the configuration file. The missing entry that triggered this error was "  + str(offender) + ".")
        sys.exit(1)

## READ IN WHETHER THE USER DOES NOT WANT US TO REWRITE THE CLUSTER CONFIGURATION FILE ##
WRITECLUSTER = 1
if "rewritecluster" in config:
    if config["rewritecluster"] == 0:
        WRITECLUSTER = 0

## READ IN WHETHER WE WOULD LIKE EUKRHYTHMIC TO CONTINUE INCOMPLETE ASSEMBLIES ##
CONTINUEFLAG = 1
if "continue" in config:
    if config["continue"] == 0:
        CONTINUEFLAG = 0

## READ IN THE AVERAGE READ LENGTH ##
AVGREADLEN = 75
if "avgreadlen" in config:
    AVGREADLEN = config["avgreadlen"]

## READ IN WHETHER WE WOULD LIKE TO REMOVE RRNAs ##
REMOVERRNA = 0
if "removerrna" in config:
    if config["removerrna"] == 1:
        REMOVERRNA = 1    

## READ IN WHETHER WE WOULD LIKE TO FILTER SEQUENCES THAT DO NOT RECRUIT READS ##
FILTERLOWREADS = 0
if "filterlowreads" in config:
    if config["filterlowreads"] == 1:
        FILTERLOWREADS = 1    
        
## READ IN WHETHER WE WOULD ALSO LIKE TO DO SALMON MAPPING ON THE CAG AND FILTER THOSE READS SEPARATELY ##
FILTERLOWREADSCAG = 0
if "filterlowreadscag" in config:
    if config["filterlowreadscag"] == 1:
        FILTERLOWREADSCAG = 1  
        
## READ IN LOCATION OF EGGNOGMAPPER DATA ##
EGGNOG_DATA_LOC = "eggnog-data"
if "eggnogDIR" in config:
    EGGNOG_DATA_LOC = config["eggnogDIR"]
    
## READ IN EUKULELE DATABASE NAME ##
EUKULELE_DATABASE = "marmmetsp"
if "eukuleleDB" in config:
    EUKULELE_DATABASE = config["eukuleleDB"]
    
## READ IN LOCATION OF PREDOWNLOADED EUKULELE DATABASE ##
EUKULELE_REFERENCE_DIR = "None"
if "eukuleleDIR" in config:
    EUKULELE_REFERENCE_DIR = config["eukuleleDIR"]
    
## READ ALL RELEVANT DATA IN FROM CONFIGURATION FILE ##
if "metaT_sample" in config:
    DATAFILE = config["metaT_sample"]
else:
    DATAFILE = os.path.join("input", "sampledata.txt")

INPUTDIR = config["inputDIR"]
OUTPUTDIR = config["outputDIR"]
SCRATCHDIR = config["scratch"]
KMERVALS = list(config['kmers'])
KMERVAL = config['defaultkmer']
MINCONTIGLENGTH = config["mincontig"]
if "spikefile" in config:
    SPIKEFILE = config['spikefile']
else:
    SPIKEFILE = ""
DROPSPIKE = config['runbbmap']
if "spiketable" in config:
    SPIKETABLE = config["spiketable"]
TRANSDECODERORFSIZE = config["transdecodercutoff"]
MINCOVERAGECLUST2 = config["secondclustercutoff"]
ADAPTER = config["adapter"]
if os.path.isfile(SPIKEFILE):
    ISFILESPIKE = 1
else:
    ISFILESPIKE = 0
ASSEMBLERS = [curr.lower() for curr in list(config['assemblers'])]
SAMPLEINFO = pd.read_csv(DATAFILE, sep = "\t")
SAMPLEINFO.SampleID = [str(curr) for curr in SAMPLEINFO.SampleID]
SAMPLEINFO.AssemblyGroup = [str(curr) for curr in SAMPLEINFO.AssemblyGroup]
ASSEMBLYDICT = dict(zip(list(SAMPLEINFO.AssemblyGroup), list(SAMPLEINFO.SampleName)))
samplenames = list(SAMPLEINFO.SampleID);
fastqnames = list(SAMPLEINFO.FastqFile);
TDEXTENSIONS = ["cds"]
if "usebedtools" in config:
    if config["usebedtools"] == 1:
        TDEXTENSIONS = ["cds","full.cds"]

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
        
logger.info("\033[1;35m Checking directory formatting...  \n")
# Check to make sure the user hasn't added trailing /.
for dir_curr in directories:
    try:
        if ((dir_curr[len(dir_curr)-1] == '/') | (dir_curr[len(dir_curr)-1] == '\\')):
            raise ValueError
    except ValueError:
        logger.info("Please do not add trailing slashes to input, output, scratch, assembled, or renamed directories.")
        sys.exit(1)

logger.info("\033[1;35m Setting appropriate cluster.yaml entries if specified...  \n")
for r in cluster.keys():
    if WRITECLUSTER == 1: # need to add additional flag to make sure only done once
        if "account" in cluster[r].keys():
            cluster[r]["account"] = DEFAULTUSER
        if "queue" in cluster[r].keys():
            cluster[r]["queue"] = DEFAULTQUEUE
        if ("nodes" in cluster[r].keys()) & (r == "trinity"):
            cluster[r]["nodes"] = MAXNODES
        if ("mem" in cluster[r].keys()) & (r == "trinity"):
            cluster[r]["mem"] = MAXMEMORY
        if ("threads" in cluster[r].keys()) & (r == "trinity"):
            cluster[r]["threads"] = MAXTHREADS
            
if WRITECLUSTER == 1:
    with open('cluster.yaml', 'w') as f:
        yaml.dump(cluster, f)

logger.info("\033[1;35m Checking that appropriate input files exist...  \n")
foundnames = []
for filename in fastqnames:
    search_dir = INPUTDIR
    if "/" in str(filename):
        innerdir = filename.split("/")[0]
        search_dir = os.path.join(search_dir,innerdir)
    foundnames.extend(os.listdir(search_dir))

inputfiles = "|".join(list(set(foundnames)))
filenames = []
singleorpaired = []
for currfile_ind in range(0, len(fastqnames)):
    currfile = str(str(fastqnames[currfile_ind]).split("/")[-1])
    occurrences = inputfiles.count(currfile)
    if occurrences > 2:
        logger.info("There are too many occurrences of fastq file " + currfile + " in input directory.")
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
assemblygroups=[str(curr) for curr in assemblygroups]
