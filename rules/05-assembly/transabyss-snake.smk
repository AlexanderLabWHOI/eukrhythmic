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
    
def get_samples_commas_TA(assemblygroup, dropspike, leftorright, commas = False):
    samplelist = list(set(SAMPLEINFO.loc[SAMPLEINFO['AssemblyGroup'] == assemblygroup]['SampleID']))
    foldername = os.path.join("intermediate-files", "01-setup",\
                          "03-alignment-spike")
    extensionname = "clean"
    if dropspike == 0:
        foldername = os.path.join("intermediate-files", "01-setup",\
                          "02-trim")
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
    
rule transabyss:
    input:
        left = lambda filename: get_samples_commas_TA(filename.assembly, DROPSPIKE, "left", commas = False),
        right = lambda filename: get_samples_commas_TA(filename.assembly, DROPSPIKE, "right", commas = False)
    output:
         os.path.join(ASSEMBLEDDIR, "05c-transabyss", "TA-transabyss",\
                              "TA-{assembly}_{k}_transabyss.fasta-final.fa")
    params:
        extra = "",
        kval = "{k}",
        fastaname = "TA-{assembly}_{k}_transabyss.fasta",
        outdir = os.path.join(ASSEMBLEDDIR, "05c-transabyss", "TA-transabyss")
    threads: 4
    log:
        err = os.path.join(OUTPUTDIR, "logs",
                            "05-assembly", "05c-transabyss", "{assembly}", "outputlog_{assembly}_{k}.err"),
        out = os.path.join(OUTPUTDIR, "logs",
                            "05-assembly", "05c-transabyss", "{assembly}", "outputlog_{assembly}_{k}.out")
    conda: os.path.join("..","..","envs","transabyss-env.yaml")
    shell:
        '''
        transabyss --pe {input.left} {input.right} --name {params.fastaname} --outdir {params.outdir} --kmer {params.kval} --threads 8 2> {log.err} 1> {log.out}
        '''
