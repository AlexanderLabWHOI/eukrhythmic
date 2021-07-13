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

rule mad_mmseqs:
    input: 
        infiles = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly", "11-SWAM", "merged.fasta")
    output:
        outfasta = os.path.join(OUTPUTDIR, "intermediate-files",\
                                "03-merge", "12-MAD", "MAD.fasta"),
        outtsv = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                              "12-MAD", "MAD.tsv")
    params:
        threads = 10,
        maxmemory = 30000, # -G o indicates local sequence identity.
        identityparam = 1.00,
        mincoverageshorter = MINCOVERAGECLUST2,
        mincoveragelong = 0.005,
        name_db = "MAD",
        name_intermed = "MAD_2",
        name_subdb = "MAD_3"
    log:
        err = os.path.join(OUTPUTDIR, "logs", "12-MAD", "MAD.err"),
        out = os.path.join(OUTPUTDIR, "logs", "12-MAD", "MAD.log")
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
