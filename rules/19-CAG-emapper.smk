configfile: "config.yaml"
  
import io
import os
import pathlib
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

download_data=False
if not os.path.isfile(os.path.join(EGGNOG_DATA_LOC,"eggnog_proteins.dmnd")):
    download_data=True

#envvars:
#    "TMPDIR"
tmpdir="tmk-smk"

rule cag_functional:
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
        tmpdir = os.path.join(tmpdir,"tmp_{assembly}_CAG"),
        eggnog_mapper_data = EGGNOG_DATA_LOC,
        download_data = download_data
    shell:
        '''
        if [ {params.download_data} == "True" ]; then
            mkdir -p {params.eggnog_mapper_data}
            download_eggnog_data.py -y --data_dir {params.eggnog_mapper_data}
        fi
        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}
        export EGGNOG_DATA_DIR={params.eggnog_mapper_data}
        emapper.py --override -i {input.assembly_file} --itype proteins -m diamond -o {params.prefix} --output_dir {params.outdir} --temp_dir {params.tmpdir}
        '''
