configfile: "config.yaml"

import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

def combineassemblers(assembly):
    return(" ".join([os.path.join(ASSEMBLEDDIR, assembly + "_" + curr + ".fasta") for curr in ASSEMBLERS]))
    
def combineassemblerslist(assembly):
    return([os.path.join(ASSEMBLEDDIR, assembly + "_" + curr + ".fasta") for curr in ASSEMBLERS])
 
rule metaquast_mad:
    input:
        outputassemblies = os.path.join(OUTPUTDIR, "intermediate-files",\
                                        "03-merge", "12-MAD", "MAD.fasta")
    output:
        report = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                               "15-MAD-quality", "quast", "report.tsv")
    params:
        assemblers = ",".join(ASSEMBLERS),
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                               "15-MAD-quality", "quast")
    log:
        err = os.path.join(OUTPUTDIR, "logs", "15-MAD-quality", "quast.err"),
        out = os.path.join(OUTPUTDIR, "logs", "15-MAD-quality", "quast.out")
    conda:
        os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        '''
        metaquast {input.outputassemblies} -o {params.outdir} --threads 8 --labels {params.assemblers} 2> {log.err} 1> {log.out}
        '''