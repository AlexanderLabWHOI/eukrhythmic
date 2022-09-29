configfile: "config.yaml"

import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *
    
def get_samples_commas_spades(assemblygroup, dropspike, leftorright, commas = False, retlen=False,retfirst=False,retlast=False):
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
    if retfirst:
        return samplelist[0]
    if retlast:
        return samplelist[-1]
    if retlen:
        return len(samplelist)
    if commas:
        return ",".join(samplelist)
    else:
        return samplelist
    
print(get_samples_commas_spades("SH402", DROPSPIKE, "left", commas = False))
       
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
                     "rna_{assembly}", "transcripts.fasta")
    params:
        extra = "",
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05d-rnaspades", "rna_{assembly}"),
        left = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, "left", commas = True),
        right = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, "right", commas = True),
        numsamps = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, "right", commas = False, retlen=True),
        left1 = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, "left", commas = False, retfirst=True),
        left2 = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, "left", commas = False, retlast=True),
        right1 = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, "right", commas = False, retfirst=True),
        right2 = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, "right", commas = False, retlast=True),
        maxmem = MAXMEMORY,
        CPUs = MAXCPUSPERTASK * MAXTASKS
    log:
        err = os.path.join(OUTPUTDIR, "logs",\
                           "05-assembly", "05d-rnaspades", "{assembly}",\
                           "outputlog_{assembly}_merge.err"),
        out = os.path.join(OUTPUTDIR, "logs",\
                            "05-assembly", "05d-rnaspades", "{assembly}",\
                            "outputlog_{assembly}_merge.out") 
    conda: os.path.join("..", "..", "envs", "02-assembly-env.yaml")
    shell:
        '''
        echo {params.left}
        if [ -f {params.outdir}/params.txt ]; then
            spades.py --continue -o {params.outdir} 2> {log.err} 1> {log.out}
        elif [ {params.numsamps} -eq 2 ]; then
            spades.py -m 150 -t 8 --rna --pe1-1 {params.left1} --pe1-2 {params.right1} --pe2-1 {params.left2} --pe2-2 {params.right2} -o {params.outdir} 2> {log.err} 1> {log.out}
        else
            spades.py -m 150 -t 8 --rna --pe1-1 {params.left} --pe1-2 {params.right} -o {params.outdir} 2> {log.err} 1> {log.out}
        fi
        '''
   
rule rnaspades_cleanup:
    input:
        spadesfile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05d-rnaspades",\
                     "rna_{assembly}", "transcripts.fasta")
    output:
        assembled = os.path.join(ASSEMBLEDDIR, "{assembly}_rnaspades.fasta")
    params:
        extra = ""
    shell:
        '''
        cp {input.spadesfile} {output.assembled}
        '''
