configfile: "config.yaml"

import io
import os
from os import listdir
from os.path import isfile, join
import sys
import pandas as pd
import numpy as np
import pathlib
from snakemake.exceptions import print_exception, WorkflowError                 

# Check to make sure the configuration file contains all the required entries
required_entries = ["metaT_sample","inputDIR","checkqual","spikefile","dropspike",\
                    "kmers","assemblers","jobname","adapter","separategroups","outputDIR","renamedDIR",\
                    "scratch"]
for r in required_entries:
    try:
        if r not in config:
            raise ValueError
    except ValueError:
        print("Check that you have specified values in the configuration file.")
        sys.exit(1)

if "metaT_sample" in config:
    DATAFILE = config["metaT_sample"]
else:
    DATAFILE = os.path.join("input", "sampledata.txt")
    
INPUTDIR = config["inputDIR"]
OUTPUTDIR = config["outputDIR"]
SCRATCHDIR = config["scratch"]
KMERVALS = list(config['kmers'])

# Check to make sure the scratch directory exists; otherwise, create it.
os.system("mkdir -p " + SCRATCHDIR)

# Check to make sure a list of integer k-mer values is supplied.
if isinstance(KMERVALS, list):
    for kmer in KMERVALS:
        try:
            if not isinstance(kmer + 1, int):
                raise ValueError
        except ValueError:
            print("Please only specify integer k-mer values.")
            sys.exit(1)
else:
    try:
        raise ValueError
    except ValueError:
        print("Please specify your k-mer values as a list. See examples for details.")
        sys.exit(1)

if 'renamedDIR' in config:
    RENAMEDDIR = os.path.join(OUTPUTDIR, config['renamedDIR'])
else:
    RENAMEDDIR = os.path.join(OUTPUTDIR, "renamed")
if 'assembledDIR' in config:
    ASSEMBLEDDIR = os.path.join(OUTPUTDIR, config['assembledDIR'])
else:
    ASSEMBLEDDIR = os.path.join(OUTPUTDIR, "assembled")
ASSEMBLERS = list(config['assemblers'])

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
directories = [ASSEMBLEDDIR,INPUTDIR,OUTPUTDIR,SCRATCHDIR,RENAMEDDIR]
for dir_curr in directories:
    try:
        if ((dir_curr[len(dir_curr)-1] == '/') | (dir_curr[len(dir_curr)-1] == '\\')):
            raise ValueError
    except ValueError:
        print("Please do not add trailing slashes to input, output, scratch, assembled, or renamed directories.")
        sys.exit(1)

SAMPLEINFO = pd.read_csv(DATAFILE, sep = "\t")
samplenames = list(SAMPLEINFO.SampleID);
fastqnames = list(SAMPLEINFO.FastqFile);

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

include: "modules/fastqc-snake"
include: "modules/bbmap-snake"
include: "modules/trimmomatic-snake"
include: "modules/fastqc-trimmed-snake"
include: "modules/trinity-snake"
include: "modules/velvet-snake"
include: "modules/megahit-snake"
include: "modules/transabyss-snake"
include: "modules/transabyss-merge-snake"
include: "modules/quast-snake"
include: "modules/cd-hit-snake"
include: "modules/manipnames-snake"
include: "modules/transdecoder-snake"
include: "modules/busco-snake"
include: "modules/salmon-snake"

ruleorder: trimmomatic > trimmomatic_SE
ruleorder: trinity > trinity_SE
ruleorder: megahit > megahit_SE
ruleorder: velvet > velvet_SE
    
rule all:
    input:
        # FASTQC OUTPUTS
        fastqc1 = expand(["{base}/qc/fastqc/{sample}_{num}_fastqc.html", "{base}/qc/fastqc/{sample}_{num}_fastqc.zip"], zip, base = OUTPUTDIR, sample = filenames, num = singleorpaired),
        # MULTIQC OUTPUTS
        multiqc1 = expand("{base}/qc/multiqc/firstqcreport/multiqc_report.html", zip, base = OUTPUTDIR),
        # BBMAP OUTPUTS
        bbmap = expand(os.path.join("{base}", "bbmap", "{sample}_{num}.clean.fastq.gz"), zip, base = OUTPUTDIR, sample = filenames, num = singleorpaired),
        # TRIMMOMATIC OUTPUTS
        trimmed = expand(["{base}/firsttrim/{sample}_1.trimmed.fastq.gz", "{base}/firsttrim/{sample}_2.trimmed.fastq.gz"], zip, base = OUTPUTDIR, sample = filenames),
        # FASTQC 2 OUTPUTS (trimmed)
        fastqc2 = expand(["{base}/qc/fastqc_trimmed/{sample}_{num}.trimmed_fastqc.html", "{base}/qc/fastqc_trimmed/{sample}_{num}.trimmed_fastqcs.zip"], zip, base = OUTPUTDIR, sample = filenames, num = singleorpaired),
        # MULTIQC 2 OUTPUTS
        multiqc2 = expand("{base}/qc/multiqc/trimmedqcreport/multiqc_report.html", zip, base = OUTPUTDIR),
        # ASSEMBLER OUTPUTS
        assemblersout = expand(os.path.join("{base}", "{assembly}_{assembler}.fasta"), base = ASSEMBLEDDIR, assembly = assemblygroups, assembler = ASSEMBLERS), 
        # QUAST OUTPUTS
        quast = expand(os.path.join("{base}", "quast", "{assembly}"), base = OUTPUTDIR, assembly = assemblygroups),
        # INDIVIDUAL CLUSTERING OUTPUTS
        clustering1 = expand(os.path.join("{base}", "cluster1", "{assembly}_{assembler}.fasta"), base = OUTPUTDIR, assembly = assemblygroups, assembler = ASSEMBLERS),
        # MERGED CLUSTERING OUTPUTS
        clustering2 = expand(os.path.join("{base}", "cluster2", "{assembly}_merged.fasta"), base = OUTPUTDIR, assembly = assemblygroups),
        # TRANSDECODER OUTPUTS - temporarily, run TransDecoder on the individual assemblies
        transdecoder_temp = expand(os.path.join("{base}", "transdecoder_indiv", "{assembly}_{assembler}.fasta.transdecoder.cds"),  base = OUTPUTDIR, assembly = assemblygroups, assembler = ASSEMBLERS),
        # TRANSDECODER OUTPUTS - second merging occurs within this step
        transdecoder = expand(os.path.join("{base}", "transdecoder", "{assembly}.fasta.transdecoder.cds"),  base = OUTPUTDIR, assembly = "merged"),
        # SALMON QUANTIFICATION OF RAW AGAINST INDIVIDUAL ASSEMBLIES/ASSEMBLERS
        salmon_indiv = expand(os.path.join("{base}", "salmon_indiv", "salmon_quant_assembly_{assembly}_{assembler}"), base = OUTPUTDIR, assembly = assemblygroups, assembler = ASSEMBLERS),
        # TRANSDECODED CLUSTERING OUTPUTS
        clustering3 = expand(os.path.join("{base}", "cluster3", "{assembly}_transdecoded.fasta"), base = OUTPUTDIR, assembly = "merged"),
        # SALMON QUANTIFICATION OF RAW AGAINST FINAL ASSEMBLY
        salmon = expand(os.path.join("{base}", "salmon", "salmon_quant_assembly_{assembly}"), base = OUTPUTDIR, assembly = "merged"),
        # QUAST QUALITY ASSESSMENT OF FINAL ASSEMBLY
        quastfinal = expand(os.path.join("{base}", "quastfinal", "{assembly}"), base = OUTPUTDIR, assembly = "merged"),
        quastmerged = expand(os.path.join("{base}", "quastmerged", "{assembly}"), base = OUTPUTDIR, assembly = assemblygroups),
        # BUSCO ASSESSMENT OF FINAL ASSEMBLY
        busco = expand(os.path.join("{base}", "busco", "{assembly}"), base = OUTPUTDIR, assembly = "merged")
