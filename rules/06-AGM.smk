configfile: "config.yaml"

import io
import os
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *
    
def combineassemblers_agm(assembly):
    return(" ".join([os.path.join(ASSEMBLEDDIR, assembly + "_" + curr + "_renamed.fasta") for curr in ASSEMBLERS]))
  
def combineassemblerslist_agm(assembly):
    return([os.path.join(ASSEMBLEDDIR, assembly + "_" + curr + "_renamed.fasta") for curr in ASSEMBLERS])
# this rule is meant to change the name of the headers for the FASTA files for each of the assemblies so that 
# we can track where each of the headers originally came from

rule rename:
    input:
        assemblyfiles = os.path.join(ASSEMBLEDDIR, "{assembly}_{assembler}.fasta")
    output:
        os.path.join(ASSEMBLEDDIR, "{assembly}_{assembler}_renamed.fasta")
    params:
        assemblername = "{assembler}",
        assembly = "{assembly}"
    shell:
        '''
        sed 's/>/>{params.assemblername}_{params.assembly}_/' {input.assemblyfiles} > {output}
        '''

# this rule merges the output of the first cd-hit-est clustering step
rule merge:
    input:
        assemblyfiles = lambda wildcards: combineassemblerslist_agm(wildcards.assembly)
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly", "06-AGM", "{assembly}_merged.fasta")
    params:
        assemblyfiles = lambda wildcards: combineassemblers_agm(wildcards.assembly)
    shell:
        '''
        if [ ! -f {output} ]; then
            cat {params.assemblyfiles} > {output}
        fi
        '''