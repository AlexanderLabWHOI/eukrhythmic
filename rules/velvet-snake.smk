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

def get_samples(assemblygroup):
    samplelist = list(SAMPLEINFO.loc[SAMPLEINFO['AssemblyGroup'] == assemblygroup]['SampleID']) 
    return samplelist

if DROPSPIKE == 0:
    LEFTFILE = lambda filename: expand(os.path.join(OUTPUTDIR, "firsttrim", "{samples}_1.trimmed.fastq.gz"), samples = get_samples(filename.assembly))
    RIGHTFILE = lambda filename: expand(os.path.join(OUTPUTDIR, "firsttrim", "{samples}_2.trimmed.fastq.gz"), samples = get_samples(filename.assembly))
else:
    LEFTFILE = lambda filename: expand(os.path.join(OUTPUTDIR, "bbmap", "{samples}_1.clean.fastq.gz"), samples = get_samples(filename.assembly))
    RIGHTFILE = lambda filename: expand(os.path.join(OUTPUTDIR, "bbmap", "{samples}_2.clean.fastq.gz"), samples = get_samples(filename.assembly))
    
rule velvet:
    input:
        r1 = LEFTFILE,
        r2 = RIGHTFILE
    output:
        velvetfile = os.path.join(OUTPUTDIR, "velvet", "{assembly}", "contigs.fa"),
        metavelvetfile = os.path.join(OUTPUTDIR, "velvet", "{assembly}", "meta-velvetg.contigs.fa")
    params:
        velvetdir = directory(os.path.join(OUTPUTDIR, "velvet", "{assembly}")),
        mincontig = MINCONTIGLENGTH,
        kmerval = statistics.median(KMERVALS)
    log:
        err = os.path.join("logs","velvet","outputlog_{assembly}_err.log"),
        out = os.path.join("logs","velvet","outputlog_{assembly}_out.log")
    conda: "../envs/velvet-env.yaml"
    shell:
        '''
        velveth {params.velvetdir} {params.kmerval} -strand_specific -shortPaired -fastq {input.r1} {input.r2} 2> {log.err} 1> {log.out}
        velvetg {params.velvetdir} -cov_cutoff auto -min_contig_lgth {params.mincontig} -exp_cov auto 2>> {log.err} 1>> {log.out}
        meta-velvetg {params.velvetdir} 2>> {log.err} 1>> {log.out}
        '''
        
rule velvet_SE:
    input:
        r1 = LEFTFILE
    output:
        velvetfile = os.path.join(OUTPUTDIR, "velvet", "{assembly}", "contigs.fa"),
        metavelvetfile = os.path.join(OUTPUTDIR, "velvet", "{assembly}", "meta-velvetg.contigs.fa")
    params:
        velvetdir = directory(os.path.join(OUTPUTDIR, "velvet", "{assembly}")),
        mincontig = MINCONTIGLENGTH,
        kmerval = statistics.median(KMERVALS)
    log:
        err = os.path.join("logs","velvet","outputlog_{assembly}_err.log"),
        out = os.path.join("logs","velvet","outputlog_{assembly}_out.log")
    conda: "../envs/velvet-env.yaml"
    shell:
        '''
        velveth {params.velvetdir} {params.kmerval} -strand_specific -short -fastq {input.r1} 2> {log.err} 1> {log.out}
        velvetg {params.velvetdir} -cov_cutoff auto -min_contig_lgth {params.mincontig} -exp_cov auto 2>> {log.err} 1>> {log.out}
        meta-velvetg {params.velvetdir} 2>> {log.err} 1>> {log.out}
        '''
        
rule velvet_cleanup:
    input:
        velvetfile = os.path.join(OUTPUTDIR, "velvet", "{assembly}", "meta-velvetg.contigs.fa")
    output:
        assembled = os.path.join(ASSEMBLEDDIR, "{assembly}_velvet.fasta"),
        scratchout = directory(os.path.join(SCRATCHDIR, "velvet", "{assembly}"))
    params:
        outdir = os.path.join(OUTPUTDIR, "velvet", "{assembly}"),
        scratch = os.path.join(SCRATCHDIR, "velvet")
    shell:
        '''
        mkdir -p {params.scratch}
        cp {input.velvetfile} {output.assembled}
        if [ {params.outdir} != {params.scratch} ]
        then
            mv {params.outdir} {params.scratch}
        fi
        '''
