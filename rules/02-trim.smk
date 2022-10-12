configfile: "config.yaml"

import io
import os
from os import listdir
import glob
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

def get_filenames_trim(sample, leftorright):
    filenames = os.listdir(INPUTDIR)
    matchednames = []
    indexsample = [ind for ind in range(0,len(samplenames)) if sample == samplenames[ind]][0]
    sample_fastq = fastqnames[indexsample]
    for fcurr in filenames: 
        if sample_fastq in fcurr:
            matchednames.append(fcurr)
    matchednames = sorted(matchednames)
    if leftorright == "left":
        return matchednames[0]
    else:
        if len(matchednames) > 1:
            return matchednames[1]   
        else:
            return ""
        
def get_filenames_trim(sample, leftorright):
    matchednames = []
    indexsample = [ind for ind in range(0,len(samplenames)) if sample.strip() == samplenames[ind]][0]
    sample_fastq = fastqnames[indexsample]
    
    search_dir = INPUTDIR
    if len(sample_fastq.split("/")) > 1:
        innerdir = sample_fastq.split("/")[0]
        search_dir = os.path.join(search_dir,innerdir)
        sample = sample_fastq.split("/")[1]
    filenames = os.listdir(search_dir)
    
    for fcurr in filenames: 
        if sample_fastq.split("/")[-1] in fcurr:
            matchednames.append(fcurr)
    matchednames = sorted(matchednames)
    if leftorright == "left":
        return os.path.join(search_dir,matchednames[0])
    else:
        if len(matchednames) > 1:
            return os.path.join(search_dir,matchednames[1]) 
        else:
            return ""

rule trimmomatic:
    input:
        r1 = lambda filename: expand(os.path.join("{sampnames}"),\
                                     sampnames = get_filenames_trim(filename.sample, "left")),
        r2 = lambda filename: expand(os.path.join("{sampnames}"),\
                                     sampnames = get_filenames_trim(filename.sample, "right")),
        adapters = ADAPTER
    output:
        p1 = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                          "02-trim",\
                          "{sample}_1.trimmed.fastq.gz"),
        p2 = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                          "02-trim",\
                          "{sample}_2.trimmed.fastq.gz"),
        s1 = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                          "02-trim",\
                          "{sample}_1.unpaired.fastq.gz"),
        s2 = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                          "02-trim",\
                          "{sample}_2.unpaired.fastq.gz")
    log:
        err = os.path.join(OUTPUTDIR, "logs",\ 
                           "01-setup", "02-trim",\ 
                           "PE_{sample}_err.log"),
        out = os.path.join(OUTPUTDIR, "logs",\
                           "01-setup", "02-trim",\
                           "PE_{sample}_out.log")
    conda: os.path.join("..", "envs", "01-setup-env.yaml")
    shell:
        '''
        trimmomatic PE {input.r1} {input.r2} {output.p1} {output.s1} {output.p2} {output.s2} ILLUMINACLIP:{input.adapters}:2:30:7 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:50 2> {log.err} 1> {log.out}
        '''
        
rule trimmomatic_SE:
    input:
        r1 = lambda filename: expand(os.path.join("{sampnames}"),\
                                     sampnames = get_filenames_trim(filename.sample, "left")),
        adapters = ADAPTER
    output:
        p1 = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                          "02-trim",\
                          "{sample}_1.trimmed.fastq.gz"),
        s1 = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                          "02-trim",\
                          "{sample}_1.unpaired.fastq.gz")
    log:
        err = os.path.join(OUTPUTDIR, "logs",\
                           "01-setup", "02-trim",\
                           "SE_{sample}_err.log"),
        out = os.path.join(OUTPUTDIR, "logs",\
                           "01-setup", "02-trim",\
                           "SE_{sample}_out.log")
    conda: os.path.join("..", "envs", "01-setup-env.yaml")
    shell:
        '''
        trimmomatic SE {input.r1} {output.p1} {output.s1} ILLUMINACLIP:{input.adapters}:2:30:7 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:50 2> {log.err} 1> {log.out}
        '''
