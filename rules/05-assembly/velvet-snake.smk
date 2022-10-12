configfile: "config.yaml"

import statistics
import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
import statistics
sys.path.insert(1, '../scripts')
from importworkspace import *

def get_samples_commas_velvet(assemblygroup, dropspike, filterrrnas, leftorright, commas = False):
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
    
rule velvet:
    input:
        r1 = lambda filename: get_samples_commas_velvet(filename.assembly, DROPSPIKE, REMOVERRNA, "left", commas = False),
        r2 = lambda filename: get_samples_commas_velvet(filename.assembly, DROPSPIKE, REMOVERRNA, "right", commas = False)
    output:
        velvetfile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05e-velvet", "{assembly}", "contigs.fa"),
        metavelvetfile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05e-velvet", "{assembly}", "meta-velvetg.contigs.fa")
    params:
        velvetdir = directory(os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05e-velvet", "{assembly}")),
        mincontig = MINCONTIGLENGTH,
        kmerval = statistics.median(KMERVALS)
    log:
        err = os.path.join(OUTPUTDIR, "logs", "05e-velvet", "{assembly}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "05e-velvet", "{assembly}.out")
    conda: os.path.join("..", "..", "envs", "02-assembly-env.yaml")
    shell:
        '''
        velveth {params.velvetdir} {params.kmerval} -strand_specific -shortPaired -fastq {input.r1} {input.r2} 2> {log.err} 1> {log.out}
        velvetg {params.velvetdir} -cov_cutoff auto -min_contig_lgth {params.mincontig} -exp_cov auto 2>> {log.err} 1>> {log.out}
        meta-velvetg {params.velvetdir} 2>> {log.err} 1>> {log.out}
        '''
        
rule velvet_SE:
    input:
        r1 = lambda filename: get_samples_commas_velvet(filename.assembly, DROPSPIKE, REMOVERRNA, "left", commas = False)
    output:
        velvetfile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05e-velvet", "{assembly}", "contigs.fa"),
        metavelvetfile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05e-velvet", "{assembly}", "meta-velvetg.contigs.fa")
    params:
        velvetdir = directory(os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05e-velvet", "{assembly}")),
        mincontig = MINCONTIGLENGTH,
        kmerval = statistics.median(KMERVALS)
    log:
        err = os.path.join(OUTPUTDIR, "logs", "05e-velvet", "{assembly}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "05e-velvet", "{assembly}.log")
    conda: os.path.join("..", "..", "envs", "02-assembly-env.yaml")
    shell:
        '''
        velveth {params.velvetdir} {params.kmerval} -strand_specific -short -fastq {input.r1} 2> {log.err} 1> {log.out}
        velvetg {params.velvetdir} -cov_cutoff auto -min_contig_lgth {params.mincontig} -exp_cov auto 2>> {log.err} 1>> {log.out}
        meta-velvetg {params.velvetdir} 2>> {log.err} 1>> {log.out}
        '''
        
rule velvet_cleanup:
    input:
        velvetfile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                                  "05-assembly", "05e-velvet", "{assembly}", "meta-velvetg.contigs.fa")
    output:
        assembled = os.path.join(ASSEMBLEDDIR, "{assembly}_velvet.fasta"),
        scratchout = directory(os.path.join(SCRATCHDIR, "05e-velvet", "{assembly}"))
    params:
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                                  "05-assembly", "05e-velvet", "{assembly}"),
        scratch = os.path.join(SCRATCHDIR, "05e-velvet")
    shell:
        '''
        mkdir -p {params.scratch}
        cp {input.velvetfile} {output.assembled}
        if [ {params.outdir} != {params.scratch} ]
        then
            mv {params.outdir} {params.scratch}
        fi
        '''
