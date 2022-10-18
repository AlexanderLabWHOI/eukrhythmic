configfile: "config.yaml"

import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *
    
def get_samples_commas(assemblygroup, dropspike, filterrrnas, leftorright, commas = False):
    samplelist = list(SAMPLEINFO.loc[SAMPLEINFO['AssemblyGroup'] == assemblygroup]['SampleID']) 
    foldername = os.path.join("intermediate-files", "01-setup",\
                          "03-alignment-spike")
    extensionname = "clean"
    if dropspike == 0:
        foldername = os.path.join("intermediate-files", "01-setup",\
                          "02-trim")
        extensionname = "trimmed"
        
    if filterrrnas == 1:
        foldername = os.path.join("intermediate-files", "01-setup",\
                          "04a-ribo")
        extensionname = "ribodetector_filt"
        
    if leftorright == "left":
        samplelist = [os.path.join(OUTPUTDIR, foldername, sample + "_1." + extensionname + ".fastq.gz") 
                      for sample in samplelist]
    else:
        samplelist = [os.path.join(OUTPUTDIR, foldername, sample + "_2." + extensionname + ".fastq.gz") 
                      for sample in samplelist]
    if commas:
        return ",".join(samplelist)
    else:
        return samplelist

MINKVAL = int(min(KMERVALS))
MAXKVAL = int(max(KMERVALS))
if MINKVAL % 2 == 0:
    MINKVAL = MINKVAL + 1
if MAXKVAL % 2 == 0:
    MAXKVAL = MAXKVAL + 1

rule megahit:
    input:
        left = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, REMOVERRNA, "left", commas = False),
        right = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, REMOVERRNA, "right", commas = False)
    output:
        megafile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                          "05-assembly", "05b-megahit", "megahit_{assembly}", "final.contigs.fa")
    params:
        megadir = directory(os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                          "05-assembly", "05b-megahit", "megahit_{assembly}")),
        left = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, REMOVERRNA, "left", commas = True),
        right = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, REMOVERRNA, "right", commas = True),
        minkval = MINKVAL,
        maxkval = MAXKVAL
    log:
        err = os.path.join(OUTPUTDIR, "logs",\
                          "05-assembly", "05b-megahit", "megahit_{assembly}", "outputlog_{assembly}.err"),
        out = os.path.join(OUTPUTDIR, "logs",\
                          "05-assembly", "05b-megahit", "megahit_{assembly}", "outputlog_{assembly}.out")
    conda: os.path.join("..", "..", "envs", "02-assembly-env.yaml")
    shell:
        '''
        if [ -d {params.megadir} ]
        then
            megahit --continue --k-min {params.minkval} --k-max {params.maxkval} -m 0.9 -t 8 -1 {params.left} -2 {params.right} -f -o {params.megadir} 2> {log.err} 1> {log.out}
        else
            megahit --k-min {params.minkval} --k-max {params.maxkval} -m 0.9 -t 8 -1 {params.left} -2 {params.right} -f -o {params.megadir} 2> {log.err} 1> {log.out}
        fi
        '''
        
rule megahit_SE:
    input:
        left = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, REMOVERRNA, "left", commas = False)
    output:
        megafile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                          "05-assembly", "05b-megahit", "megahit_{assembly}", "final.contigs.fa")
    params:
        megadir = directory(os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                          "05-assembly", "05b-megahit", "megahit_{assembly}")), 
        left = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, REMOVERRNA, "left", commas = True),
        minkval = MINKVAL,
        maxkval = MAXKVAL
    log:
        err = os.path.join(OUTPUTDIR, "logs",\
                          "05-assembly", "05b-megahit", "megahit_{assembly}", "outputlog_{assembly.err"),
        out = os.path.join(OUTPUTDIR, "logs",\
                          "05-assembly", "05b-megahit", "megahit_{assembly}", "outputlog_{assembly}.out")
    conda: os.path.join("..", "..", "envs", "02-assembly-env.yaml")
    shell:
        '''
        if [ -d {params.megadir} ]
        then
            megahit --continue --k-min {params.minkval} --k-max {params.maxkval} -m 0.9 -t 8 -r {input.left} -f -o {params.megadir} 2> {log.err} 1> {log.out}
        else
            megahit --k-min {params.minkval} --k-max {params.maxkval} -m 0.9 -t 8 -r {input.left} -f -o {params.megadir} 2> {log.err} 1> {log.out}
        fi
        '''
        
rule megahit_cleanup:
    input:
        megahitfile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                          "05-assembly", "05b-megahit", "megahit_{assembly}", "final.contigs.fa")
    output:
        assembled = os.path.join(ASSEMBLEDDIR, "{assembly}_megahit.fasta")
    params:
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                          "05-assembly", "05b-megahit", "megahit_{assembly}"),
        scratch = os.path.join(SCRATCHDIR, "02-assembly",\
                          "05-assembly", "05b-megahit", "{assembly}")
    shell:
        '''
        mkdir -p {params.scratch}
        cp {input.megahitfile} {output.assembled}
        if [ {params.outdir} != {params.scratch} ]
        then
            mv {params.outdir} {params.scratch}
        fi
        '''
