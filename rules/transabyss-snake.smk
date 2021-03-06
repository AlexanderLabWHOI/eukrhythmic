configfile: "config.yaml"

import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

def get_samples(assemblygroup):
    samplelist = list(SAMPLEINFO.loc[SAMPLEINFO['AssemblyGroup'] == assemblygroup]['SampleID']) 
    return samplelist

if DROPSPIKE == 0:
    LEFTFILE = lambda filename: expand(os.path.join(OUTPUTDIR, "firsttrim", "{samples}_1.trimmed.fastq.gz"), samples = get_samples(filename.assembly))
    RIGHTFILE = lambda filename: expand(os.path.join(OUTPUTDIR, "firsttrim", "{samples}_2.trimmed.fastq.gz"), samples = get_samples(filename.assembly))
else:
    LEFTFILE = lambda filename: expand(os.path.join(OUTPUTDIR, "bbmap", "{samples}_1.clean.fastq.gz"), samples = get_samples(filename.assembly))
    RIGHTFILE = lambda filename: expand(os.path.join(OUTPUTDIR, "bbmap", "{samples}_2.clean.fastq.gz"), samples = get_samples(filename.assembly))
    
rule transabyss:
    input:
        left = LEFTFILE,
        right = RIGHTFILE
    output:
        os.path.join(OUTPUTDIR, "transabyss_{k}_{assembly}", "{assembly}_{k}_transabyss.fasta-final.fa")
    params:
        extra = "",
        kval = "{k}",
        fastaname = "{assembly}_{k}_transabyss.fasta",
        outdir = os.path.join(OUTPUTDIR, "transabyss_{k}_{assembly}"),
        left = LEFTFILE,
        right = RIGHTFILE
    threads: 4
    log:
        err = os.path.join("logs","transabyss","outputlog_{assembly}_{k}_err.log"),
        out = os.path.join("logs","transabyss","outputlog_{assembly}_{k}_out.log")
    conda: '../envs/transabyss-env.yaml'
    shell:
        '''
        transabyss --pe {input.left} {input.right} --name {params.fastaname} --outdir {params.outdir} --kmer {params.kval} --threads 8 2> {log.err} 1> {log.out}
        '''
