configfile: "config.yaml"

import io
import os
import pathlib
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

rule link_reads:
    input:
        assembly_file = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "13-MAD-proteins-{filter_keyword}", "MAD.fasta.transdecoder.pep"),
        salmon_samples = expand(os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",
                                             "14-MAD-mapping-{filter_keyword}",
                                             "salmon_sample", "{sample}_quant",
                                             "quant.sf"), sample=filenames, filter_keyword="{filter_keyword}")
    shell:
        """
        
        """
    
rule link_function:
    input:
        hits_file = os.path.join(OUTPUTDIR, "intermediate-files",
                                 "04-compare", "17-MAD-emapper-{filter_keyword}",
                                 "MAD.emapper.hits")
    output:
    
    params:
        annotation_file = os.path.join(OUTPUTDIR, "intermediate-files",
                                 "04-compare", "17-MAD-emapper-{filter_keyword}",
                                 "MAD.emapper.annotations")