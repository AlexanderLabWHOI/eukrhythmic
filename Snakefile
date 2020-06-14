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
sys.path.insert(1, 'scripts')

# Import relevant variable from config; exit if not supplied
import importworkspace
from importworkspace import *
os.system("python " + os.path.join("scripts", "importworkspace.py"))

# Contains function to check that variables are present and formatted correctly beyond provided values in config file
import checkrequirements
from checkrequirements import *

## CHECK THAT ALL REQUIREMENTS ARE SATISFIED BY EXECUTING checkrequirements() from `scripts/checkrequirements.py` ##
checkrequirementsfct()

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
include: "modules/hardclean-snake"

ruleorder: trimmomatic > trimmomatic_SE
ruleorder: trinity > trinity_SE
ruleorder: megahit > megahit_SE
ruleorder: velvet > velvet_SE

rule all:
    input:
        # FASTQC OUTPUTS
        fastqc1 = expand([os.path.join("{base}", "qc", "fastqc", "{sample}_{num}_fastqc.html"), os.path.join("{base}", "qc", "fastqc", "{sample}_{num}_fastqc.zip")], zip, base = OUTPUTDIR, sample = filenames, num = singleorpaired),
        # MULTIQC OUTPUTS
        multiqc1 = expand(os.path.join("{base}", "qc", "multiqc", "firstqcreport", "multiqc_report.html"), zip, base = OUTPUTDIR),
        # BBMAP OUTPUTS
        bbmap = expand(os.path.join("{base}", "bbmap", "{sample}_{num}.clean.fastq.gz"), zip, base = OUTPUTDIR, sample = filenames, num = singleorpaired),
        # TRIMMOMATIC OUTPUTS
        trimmed = expand([os.path.join("{base}", "firsttrim", "{sample}_1.trimmed.fastq.gz"), os.path.join("{base}", "firsttrim", "{sample}_2.trimmed.fastq.gz")], zip, base = OUTPUTDIR, sample = filenames),
        # FASTQC 2 OUTPUTS (trimmed)
        fastqc2 = expand([os.path.join("{base}", "qc", "fastqc_trimmed", "{sample}_{num}.trimmed_fastqc.html"), os.path.join("{base}", "qc", "fastqc_trimmed", "{sample}_{num}.trimmed_fastqc.zip")], zip, base = OUTPUTDIR, sample = filenames, num = singleorpaired),
        # MULTIQC 2 OUTPUTS
        multiqc2 = expand(os.path.join("{base}", "qc", "multiqc", "trimmedqcreport", "multiqc_report.html"), zip, base = OUTPUTDIR),
        # ASSEMBLER OUTPUTS
        assemblersout = expand(os.path.join("{base}", "{assembly}_{assembler}.fasta"), base = ASSEMBLEDDIR, assembly = assemblygroups, assembler = ASSEMBLERS), 
        # QUAST OUTPUTS
        quast = expand(os.path.join("{base}", "quast", "{assembly}"), base = OUTPUTDIR, assembly = assemblygroups),
        # COMBINE QUAST OUTPUTS
        quastcombine = expand(os.path.join("{base}", "quast", "fullresults", "allresults.tsv"), base = OUTPUTDIR),
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
        # COMBINE QUAST MERGED OUTPUTS
        quastmergedcombine = expand(os.path.join("{base}", "quastmerged", "fullresults", "allresults.tsv"), base = OUTPUTDIR),
        # BUSCO ASSESSMENT OF FINAL ASSEMBLY
        busco = expand(os.path.join("{base}", "busco", "{assembly}"), base = OUTPUTDIR, assembly = "merged")
        
