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
        pep = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "{filter_workflow}",\
                           "13-MAD-proteins", "MAD.fasta.transdecoder.pep")
    output:
        eukulele_done = os.path.join(OUTPUTDIR, "intermediate-files","04-compare","{filter_workflow}",
                                     "16-MAD-taxonomy",
                                     "EUKulele_done.txt"),
        taxonomy_estimation = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",
                                     "{filter_workflow}",
                                     "16-MAD-taxonomy","taxonomy_estimation",
                                     "MAD.fasta.transdecoder-estimated-taxonomy.out")
    params:
        sampledir = os.path.join(OUTPUTDIR, "intermediate-files",
                                 "04-compare",
                                 "{filter_workflow}",
                                 "13-MAD-proteins"),
        eukulele_directory = os.path.join(OUTPUTDIR, "intermediate-files", 
                                          "04-compare",
                                          "{filter_workflow}",
                                          "16-MAD-taxonomy"),
        eukulele_reference_dir = EUKULELE_REFERENCE_DIR,
        eukulele_database = EUKULELE_DATABASE
    shell:
        """
        if [ {params.eukulele_reference_dir} != "None" ]; then
            EUKulele --mets_or_mags mets --sample_dir {params.sampledir} --p_ext ".pep" --reference_dir {params.eukulele_reference_dir} -o {params.eukulele_directory}
        else
            EUKulele --mets_or_mags mets --sample_dir {params.sampledir} --p_ext ".pep" --database {params.eukulele_database} -o {params.eukulele_directory}
        fi
        touch {output.eukulele_done}
        """
