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
        extensionname = "ribodetector_rrna_reads"
        
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
rule trinity:
    input:
        left = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE,\
                                                   REMOVERRNA, "left", commas = False),
        right = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE,\
                                                    REMOVERRNA, "right", commas = False)
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                          "05-assembly", "05a-trinity", "trinity_{assembly}.Trinity.fasta")
    params:
        extra = "",
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                          "05-assembly", "05a-trinity", "trinity_{assembly}"),
        left = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, REMOVERRNA,
                                                   "left", commas = True),
        right = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, REMOVERRNA,
                                                    "right", commas = True),
        maxmem = MAXMEMORY,
        CPUs = MAXCPUSPERTASK * MAXTASKS
    log:
        err = os.path.join(OUTPUTDIR, "logs", "02-assembly",\
                           "05-assembly", "05a-trinity", "outputlog_{assembly}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "02-assembly",\
                           "05-assembly", "05a-trinity", "outputlog_{assembly}.log")
    conda: os.path.join("..", "..", "envs", "02-assembly-env.yaml")
    shell:
        '''
        echo {params.left}
        Trinity --seqType fq --min_kmer_cov 2 --max_memory {params.maxmem}G --CPU {params.CPUs} --left {params.left} --right {params.right} --output {params.outdir} --NO_SEQTK 2> {log.err} 1> {log.out}
        '''
        
rule trinity_SE:
    input:
        single = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, REMOVERRNA, 
                                                     "left", commas = False)
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                          "05-assembly", "05a-trinity", "trinity_{assembly}", "Trinity.fasta")
    params:
        extra = "",
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                          "05-assembly", "05a-trinity", "trinity_{assembly}"),
        single = lambda filename: get_samples_commas(filename.assembly, DROPSPIKE, REMOVERRNA,
                                                     "left", commas = True),
        maxmem = MAXMEMORY,
        CPUs = MAXCPUSPERTASK * MAXTASKS
    log:
        err = os.path.join(OUTPUTDIR, "logs", "02-assembly",\
                           "05-assembly", "05a-trinity", "outputlog_{assembly}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "02-assembly",\
                           "05-assembly", "05a-trinity", "outputlog_{assembly}.log")
    conda: os.path.join("..", "envs", "02-assembly-env.yaml")
    shell:
        '''
        Trinity --seqType fq --max_memory {params.maxmem}G --CPU {params.CPUs} --max_memory 150G --bflyCalculateCPU --single {params.single} --output {params.outdir} --NO_SEQTK 2> {log.err} 1> {log.out}
        '''
   
rule trinity_cleanup:
    input:
        trinityfile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                          "05-assembly", "05a-trinity", "trinity_{assembly}.Trinity.fasta")
    output:
        assembled = os.path.join(ASSEMBLEDDIR, "{assembly}_trinity.fasta")
    params:
        extra = "",
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                          "05-assembly", "05a-trinity", "trinity_{assembly}"),
        scratch = os.path.join(SCRATCHDIR, "02-assembly", "05-assembly", "05a-trinity", "trinity_{assembly}"),
        jellyfile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                          "05-assembly", "05a-trinity", "trinity_{assembly}", "jellyfish.kmers.25.asm.fa"),
        jellyfish_final = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                                 "05-assembly", "05a-trinity-jellyfish", "jellyfish_25.fasta"),
        jellyfish_dir =  os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                                 "05-assembly", "05a-trinity-jellyfish")
    shell:
        '''
        mkdir -p {params.scratch}
        cp {input.trinityfile} {output.assembled}
        mkdir -p {params.jellyfish_dir}
        #mv {params.jellyfile} {params.jellyfish_final}
        #if [ {params.outdir} != {params.scratch} ]
        #then
        #    mv {params.outdir} {params.scratch}
        #fi
        '''
