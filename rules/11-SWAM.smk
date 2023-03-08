configfile: "config.yaml"

import io
import os
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

def getallmergedlist(filter_key):
    if filter_key != "filtered":
        mergednames = list(set([os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "07-CAG",\
                                             curr + "_merged.fasta") \
                                for curr in SAMPLEINFO.AssemblyGroup]))
    else:
        mergednames = list(set([os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "19-CAG-filtered",\
                                             curr + ".filtered.fasta") \
                                for curr in SAMPLEINFO.AssemblyGroup]))
    return(mergednames)
    
def getallmerged(filter_key):
    if filter_key != "filtered":
        mergednames = list(set([os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "07-CAG",\
                                             curr + "_merged.fasta") \
                                for curr in SAMPLEINFO.AssemblyGroup]))
    else:
        mergednames = list(set([os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "19-CAG-filtered",\
                                             curr + ".filtered.fasta") \
                                for curr in SAMPLEINFO.AssemblyGroup]))
    return(" ".join(mergednames))
    
# this rule merges all of the different assemblies of samples
rule merge_all:
    input:
        assemblyfiles = lambda filename: getallmergedlist(filename.filter_key)
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly", "11-SWAM", "{filter_key}.fasta")
    params:
        assemblyfiles = lambda filename: getallmerged(filename.filter_key)
    shell:
        '''
        echo "" > {output}
        for assemblyfile in {params.assemblyfiles}; do
            currvar=$(echo $assemblyfile | rev | cut -f1 -d\/ | rev)
            sed 's/>.*/&_'"$currvar"'/' $assemblyfile >> {output}
        done
        #cat {params.assemblyfiles} > {output}
        '''
        