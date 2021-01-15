configfile: "config.yaml"

import io
import os
import pathlib
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError

SCRATCHDIR = config["scratch"]

# Rule to remove the specified scratch directory. Invoke only after having previously invoked rule all. 
rule hardclean:
    params:
        scratch = SCRATCHDIR
    shell:
        '''
        if [ -d  {params.scratch} ] 
        then
            rm -r {params.scratch}
        fi
        '''