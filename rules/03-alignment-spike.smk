configfile: "config.yaml"

import io
import os
from os.path import isfile, join
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

rule bbmap:
    input:
        p1 = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                          "02-trim",\
                          "{sample}_1.trimmed.fastq.gz"),
        p2 = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                          "02-trim",\
                          "{sample}_2.trimmed.fastq.gz")
    output:
        r1 = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                          "03-alignment-spike", "{sample}_1.clean.fastq.gz"),
        r2 = os.path.join(OUTPUTDIR, "intermediate-files", "01-setup",\
                          "03-alignment-spike", "{sample}_2.clean.fastq.gz")
    log:
        os.path.join(OUTPUTDIR, "logs", "01-setup", "03-alignment-spike",\
                     "outputlog_{sample}_bbmap.log")
    params:
        spikefile = SPIKEFILE,
        isspike = ISFILESPIKE
    threads: 4
    conda: os.path.join("..", "envs", "01-setup-env.yaml")
    shell:
        '''
        if {params.isspike}  ; then
            bbmap.sh in1={input.p1} in2={input.p2} ref={params.spikefile} outu1={output.r1} outu2={output.r2}
        else
            cp {input.p1} {output.r1}
            cp {input.p2} {output.r2}
        fi
        '''

