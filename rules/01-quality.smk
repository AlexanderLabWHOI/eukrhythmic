configfile: "config.yaml"

import io
import os
import pathlib
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

def get_filenames(sample, leftorright, stub=False):
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
    if stub:
        search_dir = ""
    if int(leftorright.strip()) == 1:
        return os.path.join(search_dir,matchednames[0])
    else:
        return os.path.join(search_dir,matchednames[1])

def splitfilenames(inname):
    return inname.split(".")[0]

def getalloutputs():
    listqc = []
    for s in SAMPLEINFO.SampleID:
        listqc.extend([os.path.join(OUTPUTDIR, "qc", "fastqc", s + "_1_fastqc.html"),\
        os.path.join(OUTPUTDIR, "qc", "fastqc", s + "_2_fastqc.html")])
    return listqc

def getalloutputs():
    listqc = []
    for s in SAMPLEINFO.SampleID:
        listqc.extend([os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                                    "01-quality", "01a-fastqc",\
                                    s + "_1_fastqc.html"),\
                       os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                                    "01-quality",  "01a-fastqc",\
                                    s + "_2_fastqc.html")])
    return listqc
    
rule fastqc:
    input:
        (lambda filename: expand(os.path.join("{sampnames}"),\
        sampnames = get_filenames(filename.sample, filename.num)))
    output:
        html = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                            "01-quality", "01a-fastqc", "{sample}_{num}_fastqc.html"),
        zip = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                           "01-quality", "01a-fastqc", "{sample}_{num}_fastqc.zip")    
    params:
        filenum = "{num}",
        fastqdir = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                                "01-quality", "01a-fastqc"),
        outname = (lambda filename: os.path.join(OUTPUTDIR, "intermediate-files",\
                  "01-setup","01-quality",\
                  "01a-fastqc", splitfilenames(get_filenames(filename.sample, \
                  filename.num, stub=True))))
    conda: os.path.join("..", "envs", "01-setup-env.yaml")
    log: 
        err = os.path.join(OUTPUTDIR, "logs", "01-setup", "01-quality",\
                           "01a-fastqc", "{sample}_{num}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "01-setup", "01-quality",\
                           "01a-fastqc", "{sample}_{num}.out")
    shell:
        '''
        mkdir -p {params.fastqdir}
        (fastqc {input} -o {params.fastqdir} 2> {log.err} 1> {log.out} ) || (true)
        if [ {params.outname}_fastqc.html != {output.html} ]; then
            mv {params.outname}_fastqc.html {output.html} 2>/dev/null
            mv {params.outname}_fastqc.zip {output.zip} 2>/dev/null
        fi
        '''
        
rule multiqc:
    input:
        fastqcfiles = getalloutputs()
    output:
        htmlreport = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup", "01-quality",\
                     "01b-multiqc", "multiqc_report.html")
    params:
        fastqcdir = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup", "01-quality",\
                    "01a-fastqc"),
        multiqc = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup", "01-quality",\
                  "01b-multiqc")
    conda: os.path.join("..", "envs", "01-setup-env.yaml")
    shell:
        '''
        export LC_ALL=en_US.utf-8
        export LANG=en_US.utf-8
        multiqc -o {params.multiqc} -f {params.fastqcdir}
        '''
