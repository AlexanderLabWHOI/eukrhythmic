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

rule clustering_by_assembly_group_mmseqs:
    input: 
        infiles = os.path.join(OUTPUTDIR, "merged", "{assembly}_merged.fasta")
    output:
        outfasta = os.path.join(OUTPUTDIR, "cluster_{folder}", "{assembly}_merged.fasta"),
        outtsv = os.path.join(OUTPUTDIR, "cluster_{folder}", "{assembly}_merged_clustered.tsv")
    params:
        threads = 10,
        maxmemory = 15000, # -G o indicates local sequence identity.
        identityparam = 1.00,
        mincoverageshorter = MINCOVERAGECLUST2,
        mincoveragelong = 0.005,
        name_db = "byassemblygroup_{assembly}",
        name_intermed = "byassemblygroup_clustered_{assembly}",
        name_subdb = "byassemblygroup_subdb_{assembly}"
    log:
        err = os.path.join("logs","mmseq","{folder}_{assembly}_err.log"),
        out = os.path.join("logs","mmseq","{folder}_{assembly}_out.log")
    conda:
        "../envs/mmseq-env.yaml"
    shell:
        '''
        mmseqs createdb {input.infiles} {params.name_db} 
        mmseqs linclust {params.name_db} {params.name_intermed} tmp --min-seq-id {params.identityparam} --cov-mode 1 -c {params.mincoverageshorter} --remove-temp-files 2> {log.err} 1> {log.out}
        mmseqs createsubdb {params.name_intermed} {params.name_db} {params.name_subdb}
        mmseqs convert2fasta {params.name_subdb} {output.outfasta}
        mmseqs createtsv {params.name_db} {params.name_db} {params.name_intermed} {output.outtsv}
        rm -f {params.name_db}*
        rm -f {params.name_intermed}*
        rm -f {params.name_subdb}*
        '''
        
rule clustering_mega_merge_mmseqs:
    input: 
        infiles = os.path.join(OUTPUTDIR, "merged_all", "merged.fasta")
    output:
        outfasta = os.path.join(OUTPUTDIR, "cluster_{folder}", "merged_merged.fasta"),
        outtsv = os.path.join(OUTPUTDIR, "cluster_{folder}", "merged_merged.fasta_clustered.tsv")
    params:
        threads = 10,
        maxmemory = 30000, # -G o indicates local sequence identity.
        identityparam = 1.00,
        mincoverageshorter = MINCOVERAGECLUST2,
        mincoveragelong = 0.005,
        name_db = "megamerge",
        name_intermed = "megamerge_clustered",
        name_intermed_labeled = "megamerge_clustered_seq",
        name_subdb = "megamerge_subdb"
    log:
        err = os.path.join("logs","mmseq","{folder}_err.log"),
        out = os.path.join("logs","mmseq","{folder}_out.log")
    conda:
        "../envs/mmseq-env.yaml"
    shell:
        '''
        mmseqs createdb {input.infiles} {params.name_db} 
        mmseqs linclust {params.name_db} {params.name_intermed} tmp --min-seq-id {params.identityparam} --cov-mode 1 -c {params.mincoverageshorter} --remove-temp-files 2> {log.err} 1> {log.out}
        mmseqs createsubdb {params.name_intermed} {params.name_db} {params.name_subdb}
        mmseqs convert2fasta {params.name_subdb} {output.outfasta}
        mmseqs createtsv {params.name_db} {params.name_db} {params.name_intermed} {output.outtsv}
        rm -f {params.name_db}*
        rm -f {params.name_intermed}*
        rm -f {params.name_subdb}*
        '''
