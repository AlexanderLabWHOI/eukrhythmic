configfile: "config.yaml"

import io
import os
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

def getallmergedlist():
    mergednames = list(set([os.path.join(OUTPUTDIR, "03-merge", "07-CAG",\
                                         "{folder}", curr + "_merged.fasta") \
                            for curr in SAMPLEINFO.AssemblyGroup]))
    return(mergednames)
    
def getallmerged():
    mergednames = list(set([os.path.join(OUTPUTDIR, "03-merge", "07-CAG",\
                                         "{folder}", curr + "_merged.fasta") \
                            for curr in SAMPLEINFO.AssemblyGroup]))
    return(" ".join(mergednames))
    
# this rule merges all of the different assemblies of samples
rule merge_all:
    input:
        assemblyfiles = getallmergedlist()
    output:
        os.path.join(OUTPUTDIR, "02-assembly", "11-SWAM", "merged.fasta")
    params:
        assemblyfiles = getallmerged()
    shell:
        '''
        cat {params.assemblyfiles} > {output}
        '''