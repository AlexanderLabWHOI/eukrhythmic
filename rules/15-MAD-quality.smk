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
 
rule metaquast:
    input:
        outputassemblies = lambda wildcards: combineassemblerslist(wildcards.assembly) 
    output:
        directory(os.path.join(OUTPUTDIR, "04-compare", "15-MAD-quality", "{assembly}"))
    params:
        assemblers = ",".join(ASSEMBLERS),
        outputassemblies = lambda wildcards: combineassemblers(wildcards.assembly) 
    log:
        err = os.path.join(OUTPUTDIR, "logs", "15-MAD-quality", "{assembly}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "15-MAD-quality", "{assembly}.out")
    conda:
        os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        '''
        python metaquast {params.outputassemblies} -o {output} --threads 8 --labels {params.assemblers} 2> {log.err} 1> {log.out}
        '''
        
rule combinequast:
    input:
        quastdir = [os.path.join(OUTPUTDIR, "04-compare", "10-CAG-quality", curr) for curr in assemblygroups]
    output:
        os.path.join(OUTPUTDIR, "04-compare", "10-CAG-quality", "combined", "all.tsv")
    params:
        outputfile = os.path.join(OUTPUTDIR, "04-compare", "10-CAG-quality", "combined", "all.tsv"),
        assemblers = ",".join(ASSEMBLERS),
        quastdir = os.path.join(OUTPUTDIR, "04-compare", "10-CAG-quality") 
    run:
        quastfiles = os.listdir(params.quastdir)
        if len(quastfiles) > 0:
            outputassemblers = pd.DataFrame()
            for r in quastfiles:
                if os.path.isfile(os.path.join(params.quastdir, r, "report.tsv")):
                    currentreport = pd.read_csv(os.path.join(params.quastdir, r, "report.tsv"), index_col = 0, sep = "\t")
                    for c in range(0,len(currentreport.columns)):
                        currentreport = currentreport.rename(columns={currentreport.columns[c]: (currentreport.columns[c] + 
                                                                                                 "_" + str(r))})
                    outputassemblers = pd.concat([outputassemblers,currentreport], axis = 1)
        outputassemblers.to_csv(path_or_buf = params.outputfile, sep = "\t")