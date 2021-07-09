configfile: "config.yaml"

import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

KMERVALS = list(config['kmers'])
MINKVAL = int(min(KMERVALS))
MAXKVAL = int(max(KMERVALS))
PREFIXKS = " ".join([".k" + str(curr) for curr in KMERVALS])
    
rule transabyssmerge:
    input:
        files = [os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                              "05-assembly", "05c-transabyss",\
                              "transabyss_" + str(curr) + "_{assembly}",\
                              "{assembly}_" + str(curr) + \
                              "_transabyss.fasta-final.fa") for curr in KMERVALS]
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05c-transabyss",\
                     "{assembly}_transabyss.fasta")
    params:
        extra = "",
        minkval = MINKVAL,
        maxkval = MAXKVAL,
        prefixks = PREFIXKS
    log:
         err = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                                         "05-assembly", "05c-transabyss", "{assembly}",\
                                         "outputlog_{assembly}_merge.err"),
         out = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                                        "05-assembly", "05c-transabyss", "{assembly}",\
                                        "outputlog_{assembly}_merge.out") 
    conda: os.path.join("..","envs","02-assembly-env.yaml")
    shell:
        '''
        transabyss-merge {input.files} --mink {params.minkval} --maxk {params.maxkval} --prefixes {params.prefixks} --out {output} 2> {log.err} 1> {log.out}
        '''

rule transabyssmerge_cleanup:
    input:
        transabyssfile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05c-transabyss",\
                     "{assembly}_transabyss.fasta")
    output:
        assembled = os.path.join(ASSEMBLEDDIR, "{assembly}_transabyss.fasta")
    params:
        mergefile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05c-transabyss",\
                     "{assembly}_transabyss.fasta"),
        outdir = os.path.join(os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                              "05-assembly", "05c-transabyss", "transabyss_*_{assembly}"),
        scratch = os.path.join(SCRATCHDIR, "transabyss")
    shell:
        '''
        mkdir -p {params.scratch}
        cp {input.transabyssfile} {output.assembled}
        if [ {params.outdir} != {params.scratch} ]
        then
            mv {params.mergefile} {params.scratch}
            mv {params.outdir} {params.scratch}
        fi
        '''
