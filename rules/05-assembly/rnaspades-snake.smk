configfile: "config.yaml"

import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *
    
def get_samples_commas_spades(assemblygroup, dropspike, leftorright, commas = False):
    samplelist = list(SAMPLEINFO.loc[SAMPLEINFO['AssemblyGroup'] == assemblygroup]['SampleID']) 
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
       
# This module needs to grab all of the list of the individual files associated with the specified
# assembly group, after the scripts/make-assembly-file.py script builds said assembly groups 
# according to user specifications.  
rule rnaspades:
    input:
        left = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, "left", commas = False),
        right = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, "right", commas = False)
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05d-rnaspades",\
                     "{assembly}", "transcripts.fasta")
    params:
        extra = "",
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05d-rnaspades", "{assembly}"),
        left = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, "left", commas = True),
        right = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, "right", commas = True),
        maxmem = MAXMEMORY,
        CPUs = MAXCPUSPERTASK * MAXTASKS
    log:
        err = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                           "05-assembly", "05d-rnaspades", "{assembly}",\
                           "outputlog_{assembly}_merge.err"),
         out = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                            "05-assembly", "05d-rnaspades", "{assembly}",\
                            "outputlog_{assembly}_merge.out") 
    conda: os.path.join("..", "envs", "02-assembly-env.yaml")
    shell:
        '''
        echo {params.left}
        spades.py --rna --pe1-1 {params.left} --pe1-2 {params.right} -o {params.outdir} 2> {log.err} 1> {log.out}
        '''
   
rule rnaspades_cleanup:
    input:
        spadesfile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05d-rnaspades",\
                     "{assembly}", "transcripts.fasta")
    output:
        assembled = os.path.join(ASSEMBLEDDIR, "{assembly}_rnaspades.fasta")
    params:
        extra = ""
    shell:
        '''
        cp {input.spadesfile} {output.assembled}
        '''