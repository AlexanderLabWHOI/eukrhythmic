configfile: "config.yaml"

import io
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
import pathlib
from snakemake.exceptions import print_exception, WorkflowError                 
                                    
DATAFILE = config["metaT_sample"]
INPUTDIR = config["inputDIR"]
INPUTDIRs = config["inputDIRs"]
OUTPUTDIR = config["outputDIR"]
ASSEMBLYFILE = config["assembly"]
SCRATCHDIR = config["scratch"]
INPUTFILES = [[os.path.join(curr,f) for f in os.listdir(os.path.join(INPUTDIR, curr)) if isfile(join(os.path.join(INPUTDIR, curr), f))] for curr in INPUTDIRs.split(",")];
INPUTFILES = [item for sublist in INPUTFILES for item in sublist]
#print(INPUTFILES)

samplenames = list(pd.read_csv(DATAFILE, sep = "\t").SampleID);


# create a dictionary that contains a list with the relevant
# information about each sample: the barcode and the L code
def get_base_names():
    extensions = [] 
    filenames = []
    for i in INPUTFILES:
        split_i = i.split("_")
        extensions.append(split_i[len(split_i)-1]) # sequence.fastq or sequence.fq
        filenames.append("_".join(split_i[0:(len(split_i)-2)]))  # we also want to cut out the 1 or 2
        
    return filenames,extensions
    
filenames,extensions = get_base_names()

print(filenames)

if config["separategroups"] == 1:
    assemblygroups = list(set(pd.read_csv(ASSEMBLYFILE, sep = "\t").AssemblyGroup))
else:
    assemblygroups = [1] * len(INPUTFILES)

include: "modules/fastqc-snake"
include: "modules/trimmomatic-snake"
include: "modules/fastqc-trimmed-snake"
include: "modules/trinity-wrapper-snake"

rule all:
    input:
        # FASTQC OUTPUTS
        fastqc1 = expand(["{base}/qc/fastqc/{sample}_{oldext}.html", "{base}/qc/fastqc/{sample}_{oldext}.zip"], zip, base = OUTPUTDIR, sample = filenames, oldext = extensions),
        # TRIMMOMATIC OUTPUTS
        trimmed = expand(["{base}/firsttrim/{sample}_{oldext}_1.trimmed.fastq.gz", "{base}/firsttrim/{sample}_{oldext}_2.trimmed.fastq.gz"], zip, base = OUTPUTDIR, sample = filenames, oldext = extensions),
        # FASTQC 2 OUTPUTS (trimmed)
        fastqc2 = expand(["{base}/qc/fastqc_trimmed/{sample}_{oldext}.trimmed.html", "{base}/qc/fastqc_trimmed/{sample}_{oldext}.trimmed.zip"], zip, base = OUTPUTDIR, sample = filenames, oldext = extensions),
        # TRINITY OUTPUTS
        trinity = expand("{base}/trinity_results_assembly_{assembly}/Trinity.fasta", base = OUTPUTDIR, assembly = assemblygroups)
