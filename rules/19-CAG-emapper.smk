configfile: "config.yaml"
  
import io
import os
import pathlib
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

#envvars:
#    "TMPDIR"
tmpdir="tmk-smk"
rule emappercag:
    input:
        assembly_file = os.path.join(OUTPUTDIR, "intermediate-files",
                                             "04-compare",\
                           "08-CAG-proteins", "{assembly}_CAG.fasta.transdecoder.pep")
    output:
        hits_file = os.path.join(OUTPUTDIR, "intermediate-files",
                                 "04-compare", "19-CAG-emapper",
                     "{assembly}.emapper.hits")
    conda: "../envs/04-compare-env.yaml"
    params:
        prefix = "{assembly}",
        outdir = os.path.join(OUTPUTDIR, "intermediate-files",
                              "04-compare", "19-CAG-emapper"),
        tmpdir = os.path.join(tmpdir,"tmp_{assembly}_CAG")
    shell:
        '''
        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}
        export EGGNOG_DATA_DIR=/vortexfs1/omics/alexander/data/databases/eggnog-mapper-data/
        emapper.py --override -i {input.assembly_file} --itype proteins -m diamond -o {params.prefix} --output_dir {params.outdir} --temp_dir {params.tmpdir}
        '''
