configfile: "config.yaml"
  
import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

rule eukulele_cag_indiv:
    input:
        pep = os.path.join(OUTPUTDIR, "intermediate-files",
                         "04-compare",\
                         "08-CAG-proteins", "{assembly}_CAG.fasta.transdecoder.pep"),
        salmon = os.path.join(OUTPUTDIR, "intermediate-files",\
                 "04-compare",
                 "09-CAG-mapping",\
                 "salmon", "{assembly}_quant", "quant.sf")
    output:
        eukulele_done = os.path.join(OUTPUTDIR, "intermediate-files","04-compare","18-CAG-taxonomy",
                                     "EUKulele_{assembly}_done.txt")
    params:
        sampledir = os.path.join(OUTPUTDIR, "intermediate-files",
                                             "04-compare",\
                                 "08-CAG-proteins"),
        eukulele_directory = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare","18-CAG-taxonomy"),
        eukulele_reference_dir = EUKULELE_REFERENCE_DIR,
        eukulele_database = EUKULELE_DATABASE
    shell:
        """
        if [ {params.eukulele_reference_dir} != "None" ]; do
            EUKulele --mets_or_mags mets --sample_dir {params.sampledir} --p_ext ".pep" --reference_dir {params.eukulele_reference_dir} -o {params.eukulele_directory}
        else
            EUKulele --mets_or_mags mets --sample_dir {params.sampledir} --p_ext ".pep" --database {params.eukulele_database} -o {params.eukulele_directory}
        fi
        touch {output.eukulele_done}
        """