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

#include: "modules/fastqc-snake"
#include: "modules/bbmap-snake"
#include: "modules/trimmomatic-snake"
#include: "modules/fastqc-trimmed-snake"
#include: "modules/trinity-snake"
#include: "modules/velvet-snake"
#include: "modules/megahit-snake"
#include: "modules/transabyss-snake"
#include: "modules/transabyss-merge-snake"
include: "modules/quast-snake"
include: "modules/cd-hit-snake"
include: "modules/manipnames-snake"
include: "modules/transdecoder-snake"
include: "modules/busco-snake"
include: "modules/salmon-snake"
include: "modules/annotate-snake"
include: "modules/hardclean-snake"

#ruleorder: trimmomatic > trimmomatic_SE
#ruleorder: trinity > trinity_SE
#ruleorder: megahit > megahit_SE
#ruleorder: velvet > velvet_SE
ruleorder: transdecoder_indiv > transdecoder_final_proteins > transdecoder > transdecoder_by_assembly
ruleorder: salmon_indiv > salmon_clustering
ruleorder: combinequastmerge > quast_merged_transdecoded

rule all:
    input:
        # FASTQC OUTPUTS
        #fastqc1 = expand([os.path.join("{base}", "qc", "fastqc", "{sample}_{num}_fastqc.html"), 
        #                  os.path.join("{base}", "qc", "fastqc", "{sample}_{num}_fastqc.zip")], zip, 
        #                 base = OUTPUTDIR, sample = filenames, num = singleorpaired),
        # MULTIQC OUTPUTS
        #multiqc1 = expand(os.path.join("{base}", "qc", "multiqc", "firstqcreport", "multiqc_report.html"), zip, 
        #                  base = OUTPUTDIR),
        # BBMAP OUTPUTS
        #bbmap = expand(os.path.join("{base}", "bbmap", "{sample}_{num}.clean.fastq.gz"), zip, base = OUTPUTDIR, 
        #               sample = filenames, num = singleorpaired),
        # TRIMMOMATIC OUTPUTS
        #trimmed = expand([os.path.join("{base}", "firsttrim", "{sample}_1.trimmed.fastq.gz"), 
        #                  os.path.join("{base}", "firsttrim", "{sample}_2.trimmed.fastq.gz")], zip, 
        #                 base = OUTPUTDIR, sample = filenames),
        # FASTQC 2 OUTPUTS (trimmed)
        #fastqc2 = expand([os.path.join("{base}", "qc", "fastqc_trimmed", "{sample}_{num}.trimmed_fastqc.html"), 
        #                  os.path.join("{base}", "qc", "fastqc_trimmed", "{sample}_{num}.trimmed_fastqc.zip")], 
        #                 zip, base = OUTPUTDIR, sample = filenames, num = singleorpaired),
        # MULTIQC 2 OUTPUTS
        #multiqc2 = expand(os.path.join("{base}", "qc", "multiqc", "trimmedqcreport", "multiqc_report.html"), zip, 
        #                  base = OUTPUTDIR),
        # ASSEMBLER OUTPUTS
        #assemblersout = expand(os.path.join("{base}", "{assembly}_{assembler}.fasta"), 
        #                       base = ASSEMBLEDDIR, assembly = assemblygroups, assembler = ASSEMBLERS), 
        # QUAST OUTPUTS
        quast = expand(os.path.join("{base}", "quast", "{assembly}"), base = OUTPUTDIR, assembly = assemblygroups),
        # COMBINE QUAST OUTPUTS
        quastcombine = expand(os.path.join("{base}", "quast", "fullresults", "allresults.tsv"), base = OUTPUTDIR),
        # TRANSDECODER OUTPUTS - optionally, run TransDecoder on the individual assemblies
        transdecoder_indiv = expand(os.path.join("{base}", "transdecoder_indiv", 
                                                 "{assembly}_{assembler}.fasta.transdecoder.cds"),  
                                    base = OUTPUTDIR, assembly = assemblygroups, assembler = ASSEMBLERS),
        # INDIVIDUAL CLUSTERING OUTPUTS
        clustering_each_assembler = expand(os.path.join("{base}", "cluster1", "{assembly}_{assembler}.fasta"), 
                                           base = OUTPUTDIR, assembly = assemblygroups, assembler = ASSEMBLERS),
        # TRANSDECODER OUTPUTS - second merging ("mega-merge") occurs within this step
        #transdecoder = expand(os.path.join("{base}", "transdecoder_{folder}", "{assembly}.fasta.transdecoder.cds"), 
        #                      base = OUTPUTDIR, folder = "mega_merge", assembly = "merged"),
        # MERGED CLUSTERING OUTPUTS - cluster on merged samples by assembly group and mega
        clustering_by_assembly_group = expand(os.path.join("{base}", "cluster_{folder}", "{assembly}_merged.fasta"), 
                                              base = OUTPUTDIR, folder = "by_assembly_group", assembly = assemblygroups),
        clustering_mega_merge = expand(os.path.join("{base}", "cluster_{folder}", "{assembly}_merged.fasta"), 
                                   base = OUTPUTDIR, folder = "mega_merge", assembly = "merged"),
        # SALMON QUANTIFICATION OF RAW AGAINST INDIVIDUAL ASSEMBLIES/ASSEMBLERS
        salmon_indiv = expand(os.path.join("{base}", "indiv_salmon", "salmon_quant_assembly_{assembly}_{assembler}"), 
                              base = OUTPUTDIR, assembly = assemblygroups, assembler = ASSEMBLERS),
        salmon_merged = expand(os.path.join("{base}", "merged_salmon", "salmon_quant_assembly_{assembly}"), 
                               base = OUTPUTDIR, assembly = assemblygroups),
        # SALMON QUANTIFICATION OF RAW AGAINST MERGED BY ASSEMBLY GROUP
        salmon_by_assembly = expand(os.path.join("{base}", "salmon_{folder}", "salmon_quant_assembly_{assembly}"), 
                                    base = OUTPUTDIR, folder = "by_assembly_group", assembly = assemblygroups),
        # SALMON QUANTIFICATION OF RAW AGAINST MEGA-MERGED ASSEMBLY
        salmon_mega_merge = expand(os.path.join("{base}", "salmon_{folder}", "salmon_quant_assembly_{assembly}"), 
                                   base = OUTPUTDIR, folder = "mega_merge", assembly = "merged"),
        # TRANSDECODER ON FINAL AND MEGA-MERGED ASSEMBLY
        transdecoder_mega_merge = expand(os.path.join("{base}", "transdecoder_{folder}_finalproteins", 
                                                    "{assembly}.fasta.transdecoder.cds"), 
                                       base = OUTPUTDIR, folder = "mega_merge", assembly = "merged"),
        trandecoder_by_assembly_group = expand(os.path.join("{base}", "transdecoder_{folder}_finalproteins", 
                                                    "{assembly}.fasta.transdecoder.cds"), 
                                       base = OUTPUTDIR, folder = "by_assembly_group", assembly = assemblygroups),
        # QUAST QUALITY ASSESSMENT OF FINAL ASSEMBLY
        quastfinal = expand(os.path.join("{base}", "quast_{folder}", "{assembly}"), base = OUTPUTDIR, 
                            folder = "mega_merge", assembly = "merged"),
        quast_merged_mega_merge = expand(os.path.join("{base}", "quast_{folder}", "{assembly}"), base = OUTPUTDIR, 
                                         folder = "by_assembly_group", assembly = assemblygroups),
        # COMBINE QUAST MERGED OUTPUTS FOR BY ASSEMBLY GROUP
        quastmergedcombine = expand(os.path.join("{base}", "quast_{folder}", "fullresults", "allresults.tsv"), 
                                    base = OUTPUTDIR, folder = "by_assembly_group"),
        # BUSCO ASSESSMENT OF FINAL ASSEMBLY
        busco = expand(os.path.join("{base}", "busco", "{database}", "{folder}", "{assembly}"), base = OUTPUTDIR, 
                       database = "eukaryota", folder = "mega_merge", assembly = "merged"), # eukaryota, bacteria
        # HMMER ALIGNMENT OF FINAL ASSEMBLY BEFORE MEGA-MERGE
        hmmer = expand(os.path.join("{base}", "pfam", "{folder}", "{assembly}.tblout"), base = OUTPUTDIR, 
                       folder = "by_assembly_group", assembly = assemblygroups),
        # DIAMOND ALIGNMENT AND KEGG ANNOTATION 
        kegg = expand(os.path.join("{base}", "kegg", "{folder}", "{assembly}_kegg.csv"), base = OUTPUTDIR, 
                      folder = "by_assembly_group", assembly = assemblygroups)
        
