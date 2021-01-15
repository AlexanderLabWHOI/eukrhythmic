configfile: "config.yaml"

import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

# --auto-lineage-euk is better if you download SEPP
rule busco:
    input:
        outputassemblies = os.path.join(OUTPUTDIR, "transdecoder_{folder}_finalproteins",\
            "{assembly}.fasta.transdecoder.pep")
    output:
        directory(os.path.join("{base}", "busco", "{database}", "{folder}", "{assembly}"))
    params:
        assembly = "{assembly}"
    conda:
        "../envs/busco-env.yaml"
    shell:
        """
        #busco --list-datasets
        #-l eukaryota_odb10
        #cd ./scripts/busco
        #python3 setup.py install --user
        #cd ../..
        #python ./scripts/busco/src/busco/run_BUSCO.py -i {input.outputassemblies} --auto-lineage-euk -o {output} -c 8 -m transcriptome --long
        busco -f -i {input.outputassemblies} -l ./input/eukaryota_odb10 -o {params.assembly}_busco -c 8 -m transcriptome --long
        mv run_{params.assembly}_busco {output}
        """
