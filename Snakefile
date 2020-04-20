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
OUTPUTDIR = config["outputDIR"]
SCRATCHDIR = config["scratch"]
KMERVALS = list(config['kmers'])
ASSEMBLEDDIR = os.path.join(OUTPUTDIR, config['assembledDIR'])
ASSEMBLERS = list(config['assemblers'])
print(KMERVALS)

SAMPLEINFO = pd.read_csv(DATAFILE, sep = "\t")
samplenames = list(SAMPLEINFO.SampleID);
fastqnames = list(SAMPLEINFO.FastqFile);
print(fastqnames)
for currfile in fastqnames:
    if isfile(os.path.join(INPUTDIR, currfile + "_R1_001.fastq.gz")):
        print("yo")
    else:
        print(os.path.join(INPUTDIR, currfile + "_R1_001.fastq.gz"))
filenames = [currfile for currfile in fastqnames if isfile(os.path.join(INPUTDIR, currfile + "_R1_001.fastq.gz"))]

print(filenames)

if config["separategroups"] == 1:
    assemblygroups = list(set(SAMPLEINFO.AssemblyGroup))
else:
    assemblygroups = [1] * len(INPUTFILES)
    
print(assemblygroups)

def combineassemblers(assembly):
    return(" ".join([os.path.join(ASSEMBLEDDIR, assembly + "_" + curr + ".fasta") for curr in ASSEMBLERS]))
    
print("Combined assembly for assembly = HN004: ")
print(combineassemblers("HN004"))

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

rule all:
    input:
        # FASTQC OUTPUTS
        fastqc1 = expand(["{base}/qc/fastqc/{sample}_{num}.html", "{base}/qc/fastqc/{sample}_{num}.zip"], zip, base = OUTPUTDIR, sample = filenames, num = [1,2]),
        # BBMAP OUTPUTS
        bbmap = expand(os.path.join("{base}", "bbmap", "{sample}_{num}.clean.fastq.gz"), zip, base = OUTPUTDIR, sample = filenames, num = [1,2]),
        # TRIMMOMATIC OUTPUTS
        trimmed = expand(["{base}/firsttrim/{sample}_1.trimmed.fastq.gz", "{base}/firsttrim/{sample}_2.trimmed.fastq.gz"], zip, base = OUTPUTDIR, sample = filenames),
        # FASTQC 2 OUTPUTS (trimmed)
        fastqc2 = expand(["{base}/qc/fastqc_trimmed/{sample}_{num}.trimmed.html", "{base}/qc/fastqc_trimmed/{sample}_{num}.trimmed.zip"], zip, base = OUTPUTDIR, sample = filenames, num = [1,2]),
        # ASSEMBLER OUTPUTS
        assemblersout = expand(os.path.join("{base}", "{assembly}_{assembler}.fasta"), base = ASSEMBLEDDIR, assembly = assemblygroups, assembler = ASSEMBLERS), 
        # ~ all of this is now obsolete ~
        # TRINITY OUTPUTS
        #trinity = expand(os.path.join("{base}", "trinity_results_assembly_{assembly}", "Trinity.fasta"), base = OUTPUTDIR, assembly = assemblygroups),
        # TRANSABYSS OUTPUTS
        #transabyss = expand(os.path.join("{base}", "transabyss_{k}_{assembly}", "{assembly}_{k}_transabyss.fasta-final.fa"), zip, base = OUTPUTDIR, assembly = assemblygroups, k = KMERVALS),
        # TRANSABYSS-MERGE OUTPUTS
        # transabyssmerged = expand(os.path.join("{base}", "transabyss", "{assembly}_transabyss.fasta"), zip, base = OUTPUTDIR, assembly = assemblygroups),
        # VELVET OUTPUTS
        # velvet = expand(os.path.join("{base}", "velvet", "{assembly}"), base = OUTPUTDIR, assembly = assemblygroups),
        # MEGAHIT OUTPUTS
        # megahit = expand(os.path.join("{base}", "megahit", "{assembly}"), base = OUTPUTDIR, assembly = assemblygroups),
        # QUAST OUTPUTS
        quast = expand(os.path.join("{base}", "quast", "{assembly}"), base = OUTPUTDIR, assembly = assemblygroups),
        # INDIVIDUAL CLUSTERING OUTPUTS
        clustering1 = expand(os.path.join("{base}", "cluster1", "{assembly}_{assembler}.fasta"), zip, base = OUTPUTDIR, assembly = assemblygroups, assembler = ASSEMBLERS),
        # TRANSDECODER OUTPUTS - merging occurs within this step
        transdecoder = expand(os.path.join("{base}", "transdecoder", "{assembly}.fasta.transdecoder.cds"), zip, base = OUTPUTDIR, assembly = assemblygroups),
        # MERGED CLUSTERING OUTPUTS
        clustering2 = expand(os.path.join("{base}", "cluster2", "{assembly}_merged.fasta"), zip, base = OUTPUTDIR, assembly = assemblygroups),
        # TRANSDECODED CLUSTERING OUTPUTS
        clustering3 = expand(os.path.join("{base}", "cluster3", "{assembly}_transdecoded.fasta"), zip, base = OUTPUTDIR, assembly = assemblygroups)