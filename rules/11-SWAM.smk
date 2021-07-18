configfile: "config.yaml"

import io
import os
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

def getallmergedlist():
    mergednames = list(set([os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "07-CAG",\
                                         curr + "_merged.fasta") \
                            for curr in SAMPLEINFO.AssemblyGroup]))
    return(mergednames)
    
def getallmerged():
    mergednames = list(set([os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "07-CAG",\
                                         curr + "_merged.fasta") \
                            for curr in SAMPLEINFO.AssemblyGroup]))
    return(" ".join(mergednames))
    
# this rule merges all of the different assemblies of samples
rule merge_all:
    input:
        assemblyfiles = getallmergedlist()
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly", "11-SWAM", "merged.fasta")
    params:
        assemblyfiles = getallmerged()
    shell:
        '''
        echo "" > {output}
        for assemblyfile in {params.assemblyfiles}; do
            currvar=$(echo $assemblyfile | cut -f6 -d\/)
            sed 's/>.*/&_'"$currvar"'/' $assemblyfile >> {output}
        done
        #cat {params.assemblyfiles} > {output}
        '''