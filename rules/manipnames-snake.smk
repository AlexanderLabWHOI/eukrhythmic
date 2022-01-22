configfile: "config.yaml"

import io
import os
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

def getallmergedlist():
    mergednames = list(set([os.path.join(OUTPUTDIR, "merged", curr + "_merged.fasta") for curr in SAMPLEINFO.AssemblyGroup]))
    return(mergednames)
    
def getallmerged():
    mergednames = list(set([os.path.join(OUTPUTDIR, "merged", curr + "_merged.fasta") for curr in SAMPLEINFO.AssemblyGroup]))
    return(" ".join(mergednames))
    
def combineassemblers(assembly):
    return(" ".join([os.path.join(OUTPUTDIR, "assembled", assembly + "_" + curr + ".fasta") for curr in ASSEMBLERS]))
  
def combineassemblerslist(assembly):
    return([os.path.join(OUTPUTDIR, "assembled", assembly + "_" + curr + ".fasta") for curr in ASSEMBLERS])

# this rule is meant to change the name of the headers for the FASTA files for each of the assemblies so that 
# we can track where each of the headers originally came from

rule rename:
    input:
        assemblyfiles = os.path.join(ASSEMBLEDDIR, "{assembly}_{assembler}.fasta")
    output:
        os.path.join(RENAMEDDIR, "{assembly}_{assembler}.fasta")
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
        assemblyfiles = lambda wildcards: combineassemblerslist(wildcards.assembly)
    output:
        os.path.join(OUTPUTDIR, "merged", "{assembly}_merged.fasta")
    params:
        assemblyfiles = lambda wildcards: combineassemblers(wildcards.assembly)
    shell:
        '''
        if [ ! -f {output} ]; then
            cat {params.assemblyfiles} > {output}
        fi
        '''
        
# this rule merges all of the different assemblies of samples
rule merge_all:
    input:
        assemblyfiles = getallmergedlist()
    output:
        os.path.join(OUTPUTDIR, "merged_all", "merged.fasta")
    params:
        assemblyfiles = getallmerged()
    shell:
        '''
        cat {params.assemblyfiles} > {output}
        '''
