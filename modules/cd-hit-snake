configfile: "config.yaml"

import io
import os
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

def get_samples(assemblygroup):
    samplelist = list(ASSEMBLYFILE.loc[ASSEMBLYFILE['AssemblyGroup'] == assemblygroup]['SampleID']) 
    return samplelist

rule cdhit_individualclusters:
    input: 
        infiles = os.path.join(OUTPUTDIR, "renamed", "{assembly}_{assembler}.fasta")
    output:
        os.path.join(OUTPUTDIR, "cluster1", "{assembly}_{assembler}.fasta")
    params:
        threads = 10,
        maxmemory = 2549, # -G o indicates local sequence identity.
        identityparam = 1.00,
        mincoverageshorter = 1.00,
        mincoveragelong = 0.005
    log:
        err = os.path.join("logs","cdhit","individual_{assembly}_{assembler}_err.log"),
        out = os.path.join("logs","cdhit","individual_{assembly}_{assembler}_out.log")
    conda:
        "../envs/cd-hit-env.yaml"    
    shell:
        '''
        cd-hit-est -d 0 -i {input.infiles} -o {output} -T {params.threads} -M {params.maxmemory} -G o -c {params.identityparam} -aS {params.mincoverageshorter} -aL {params.mincoveragelong} 2> {log.err} 1> {log.out}
        '''
    
rule clustering_by_assembly_group:
    input: 
        infiles = os.path.join(OUTPUTDIR, "merged", "{assembly}_merged.fasta")
    output:
        os.path.join(OUTPUTDIR, "cluster_{folder}", "{assembly}_merged.fasta")
    params:
        threads = 10,
        maxmemory = 15000, # -G o indicates local sequence identity.
        identityparam = 1.00,
        mincoverageshorter = MINCOVERAGECLUST2,
        mincoveragelong = 0.005
    log:
        err = os.path.join("logs","cdhit","{folder}_{assembly}_err.log"),
        out = os.path.join("logs","cdhit","{folder}_{assembly}_out.log")
    conda:
        "../envs/cd-hit-env.yaml"
    shell:
        '''
        cd-hit-est -d 0 -i {input.infiles} -o {output} -T {params.threads} -M {params.maxmemory} -G o -c {params.identityparam} -aS {params.mincoverageshorter} -aL {params.mincoveragelong} 2> {log.err} 1> {log.out}
        '''
        
rule clustering_mega_merge:
    input: 
        infiles = os.path.join(OUTPUTDIR, "merged_all", "merged.fasta")
    output:
        os.path.join(OUTPUTDIR, "cluster_{folder}", "merged_merged.fasta")
    params:
        threads = 10,
        maxmemory = 30000, # -G o indicates local sequence identity.
        identityparam = 1.00,
        mincoverageshorter = MINCOVERAGECLUST2,
        mincoveragelong = 0.005
    log:
        err = os.path.join("logs","cdhit","{folder}_err.log"),
        out = os.path.join("logs","cdhit","{folder}_out.log")
    conda:
        "../envs/cd-hit-env.yaml"
    shell:
        '''
        cd-hit-est -d 0 -i {input.infiles} -o {output} -T {params.threads} -M {params.maxmemory} -G o -c {params.identityparam} -aS {params.mincoverageshorter} -aL {params.mincoveragelong} 2> {log.err} 1> {log.out}
        '''
