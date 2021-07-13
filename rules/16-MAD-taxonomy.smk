configfile: "config.yaml"

import io
import os
import pathlib
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

rule eukulele_run:
    input:
        pep = (lambda filename: expand(os.path.join(OUTPUTDIR, "transdecoder_{folder}_finalproteins", 
                                                    "{assembly}.fasta.transdecoder.pep"), 
               assembly = assemblygroups)),
        salmon = (lambda filename: expand(os.path.join(OUTPUTDIR, "salmon_mega_merge",
                                                       "salmon_quant_assembly_merged", "quant.sf"), 
               assembly = assemblygroups, folder = filename.folder))
    output:
        eukulele_directory = directory(os.path.join(OUTPUTDIR, "eukulele_{folder}")),
        tax_est = os.path.join(OUTPUTDIR, "eukulele_{folder}", "taxonomy_estimation", )
    params:
        sampledir = os.path.join(OUTPUTDIR, "transdecoder_{folder}_finalproteins"),
        salmondir = (lambda filename: expand(os.path.join(OUTPUTDIR, "salmon_mega_merge")))
    conda:
        os.path.join("..","envs","eukulele-env.yaml")
    shell:
        """
        EUKulele --mets_or_mags mets --sample_dir {params.sampledir} --p_ext ".pep" --n_ext ".cds" --database mmetsp --use_salmon_counts --salmon_dir {params.salmondir}
        """
