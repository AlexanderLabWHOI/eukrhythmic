configfile: "config.yaml"

import io
import os
import pandas as pd
import pathlib
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

def getalloutputs():
    listqc = []
    for s in SAMPLEINFO.SampleID:
        listqc.extend([os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                                    "04-quality",\
                                    s + "_1.trimmed_fastqc.html"),\
                       os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                                    "04-quality",\
                                    s + "_2.trimmed_fastqc.html")])
    return listqc

rule fastqc_trimmed:
    input:
        os.path.join(OUTPUTDIR, "intermediate-files", "01-setup", "02-trim",\
                     "{sample}_{num}.trimmed.fastq.gz")
    output:
        html = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                            "04-quality",\
                            "{sample}_{num}.trimmed_fastqc.html"),
        zip = os.path.join(OUTPUTDIR, OUTPUTDIR, "intermediate-files", "01-setup",\
                           "04-quality", "{sample}_{num}.trimmed_fastqc.zip")
    params: 
        fastqdir = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                          "04-quality"),
        outname = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                               "04-quality", "{sample}_{num}")
    conda: os.path.join("..", "envs", "01-setup-env.yaml")
    log: 
        err = os.path.join(OUTPUTDIR, "logs", OUTPUTDIR, "intermediate-files",\
                           "01-setup", "04-quality",\
                           "{sample}_{num}_trimmed_err.log"),
        out = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup", "04-quality",\
                           "{sample}_{num}_trimmed_out.log")
    shell:
        '''
        mkdir -p {params.fastqdir}
        fastqc {input} -o {params.fastqdir} 2> {log.err} 1> {log.out}
        '''
        
rule multiqc_trimmed:
    input:
        fastqcfiles = getalloutputs()
    output:
        htmlreport = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                                  "04-quality", "multiqc",\
                                  "multiqc_report.html")
    params:
        fastqcdir = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                                  "04-quality"),
        multiqc = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                                  "04-quality", "multiqc")
    conda: os.path.join("..", "envs", "01-setup-env.yaml")
    shell:
        '''
        export LC_ALL=en_US.utf-8
        export LANG=en_US.utf-8
        multiqc -o {params.multiqc} -f {params.fastqcdir}
        '''
