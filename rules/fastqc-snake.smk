configfile: "config.yaml"

import io
import os
import pathlib
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

def get_filenames(sample, leftorright):
    filenames = os.listdir(INPUTDIR)
    matchednames = []
    indexsample = [ind for ind in range(0,len(samplenames)) if sample.strip() == samplenames[ind]][0]
    sample_fastq = fastqnames[indexsample]
    for fcurr in filenames: 
        if sample_fastq in fcurr:
            matchednames.append(fcurr)
    matchednames = sorted(matchednames)
    if int(leftorright.strip()) == 1:
        return matchednames[0]
    else:
        return matchednames[1]

def splitfilenames(inname):
    return inname.split(".")[0]

def getalloutputs():
    listqc = []
    for s in SAMPLEINFO.SampleID:
        listqc.extend([os.path.join(OUTPUTDIR, "qc", "fastqc", s + "_1_fastqc.html"), os.path.join(OUTPUTDIR, "qc", "fastqc", s + "_2_fastqc.html")])
    return listqc
    
rule fastqc:
    input:
        (lambda filename: expand(os.path.join(INPUTDIR, "{sampnames}"), sampnames = get_filenames(filename.sample, filename.num)))
    output:
        html = os.path.join(OUTPUTDIR, "qc", "fastqc", "{sample}_{num}_fastqc.html"),
        zip = os.path.join(OUTPUTDIR, "qc", "fastqc", "{sample}_{num}_fastqc.zip")    
    params:
        filenum = "{num}",
        fastqdir = os.path.join(OUTPUTDIR, "qc", "fastqc"),
        outname = (lambda filename: os.path.join(OUTPUTDIR, "qc", "fastqc", splitfilenames(get_filenames(filename.sample, filename.num))))
    conda: "../envs/fastqc-env.yaml"
    log: 
        err = os.path.join("logs","fastqc","{sample}_{num}_err.log"),
        out = os.path.join("logs","fastqc","{sample}_{num}_out.log")
    shell:
        '''
        echo hello{params.filenum}hello
        if [ {params.filenum} == 1 ]; then
            echo "hello"
        fi
        mkdir -p {params.fastqdir}
        fastqc {input} -o {params.fastqdir} 2> {log.err} 1> {log.out}
        #mv {params.outname}_fastqc.html {output.html} 2>/dev/null
        #mv {params.outname}_fastqc.zip {output.zip} 2>/dev/null
        '''
        
rule multiqc:
    input:
        fastqcfiles = getalloutputs()
    output:
        htmlreport = os.path.join(OUTPUTDIR, "qc", "multiqc", "firstqcreport", "multiqc_report.html")
    params:
        fastqcdir = os.path.join(OUTPUTDIR, "qc", "fastqc"),
        multiqc = os.path.join(OUTPUTDIR, "qc", "multiqc", "firstqcreport")
    conda: '../envs/fastqc-env.yaml'
    shell:
        '''
        export LC_ALL=en_US.utf-8
        export LANG=en_US.utf-8
        multiqc -o {params.multiqc} -f {params.fastqcdir}
        '''
