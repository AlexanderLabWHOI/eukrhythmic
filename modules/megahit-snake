configfile: "config.yaml"

import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *
    
def get_samples_commas(assemblygroup, dropspike, leftorright, commas = False):
    samplelist = list(SAMPLEINFO.loc[SAMPLEINFO['AssemblyGroup'] == assemblygroup]['SampleID']) 
    foldername = "bbmap"
    extensionname = "clean"
    if dropspike == 0:
        foldername = "firsttrim"
        extensionname = "trimmed"
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
        left = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, "left", commas = False),
        right = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, "right", commas = False)
    output:
        megafile = os.path.join(OUTPUTDIR, "megahit", "{assembly}", "final.contigs.fa")
    params:
        megadir = directory(os.path.join(OUTPUTDIR, "megahit", "{assembly}")),
        left = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, "left", commas = True),
        right = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, "right", commas = True),
        minkval = MINKVAL,
        maxkval = MAXKVAL
    log:
        err = os.path.join("logs","megahit","outputlog_{assembly}_err.log"),
        out = os.path.join("logs","megahit","outputlog_{assembly}_out.log")
    conda: "../envs/megahit-env.yaml"
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
        left = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, "left", commas = False)
    output:
        megafile = os.path.join(OUTPUTDIR, "megahit", "{assembly}", "final.contigs.fa")
    params:
        megadir = directory(os.path.join(OUTPUTDIR, "megahit", "{assembly}")), 
        left = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, "left", commas = True),
        minkval = MINKVAL,
        maxkval = MAXKVAL
    log:
        err = os.path.join("logs","megahit","outputlog_{assembly}_err.log"),
        out = os.path.join("logs","megahit","outputlog_{assembly}_out.log")
    conda: "../envs/megahit-env.yaml"
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
        megahitfile = os.path.join(OUTPUTDIR, "megahit", "{assembly}", "final.contigs.fa")
    output:
        assembled = os.path.join(ASSEMBLEDDIR, "{assembly}_megahit.fasta"),
        scratchout = directory(os.path.join(SCRATCHDIR, "megahit", "{assembly}"))
    params:
        outdir = os.path.join(OUTPUTDIR, "megahit", "{assembly}"),
        scratch = os.path.join(SCRATCHDIR, "megahit")
    shell:
        '''
        mkdir -p {params.scratch}
        cp {input.megahitfile} {output.assembled}
        if [ {params.outdir} != {params.scratch} ]
        then
            mv {params.outdir} {params.scratch}
        fi
        '''
