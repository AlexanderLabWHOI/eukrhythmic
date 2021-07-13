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
        
rule clustering_agm:
    input:
        infiles = os.path.join(OUTPUTDIR, "intermediate-files",\
                               "02-assembly", "06-AGM", "{assembly}_merged.fasta")
    output:
        outfasta = os.path.join(OUTPUTDIR, "intermediate-files",\
                                "03-merge", "07-CAG",\
                                "{assembly}_merged.fasta"),
        outtsv = os.path.join(OUTPUTDIR, "intermediate-files",\
                              "03-merge", "07-CAG",\
                              "{assembly}_merged.tsv")
    params:
        threads = 10,
        maxmemory = 15000, # -G o indicates local sequence identity.
        identityparam = 1.00,
        mincoverageshorter = MINCOVERAGECLUST2,
        mincoveragelong = 0.005,
        name_db = "CAG_{assembly}",
        name_intermed = "CAG_2_{assembly}",
        name_subdb = "CAG_3_{assembly}"
    log:
        err = os.path.join(OUTPUTDIR, "logs", "07-CAG", "{assembly}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "07-CAG", "{assembly}.log")
    conda:
        os.path.join("..", "envs", "03-merge-env.yaml")
    shell:
        '''
        mmseqs createdb {input.infiles} {params.name_db} 
        mmseqs linclust {params.name_db} {params.name_intermed} tmp --min-seq-id {params.identityparam} --cov-mode 1 -c {params.mincoverageshorter} --remove-tmp-files 2> {log.err} 1> {log.out}
        mmseqs createsubdb {params.name_intermed} {params.name_db} {params.name_subdb}
        mmseqs convert2fasta {params.name_subdb} {output.outfasta}
        mmseqs createtsv {params.name_db} {params.name_db} {params.name_intermed} {output.outtsv}
        rm -f {params.name_db}*
        rm -f {params.name_intermed}*
        rm -f {params.name_subdb}*
        '''