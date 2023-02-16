configfile: "config.yaml"

import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *
    
def get_samples_multiflag_spades(assemblygroup, dropspike, filterrrnas, leftorright, flag = "--pe1"):
    samplelist = list(SAMPLEINFO.loc[SAMPLEINFO['AssemblyGroup'] == assemblygroup]['SampleID']) 
    foldername = os.path.join("intermediate-files", "01-setup",\
                          "03-alignment-spike")
    extensionname = "clean"
    if dropspike == 0:
        foldername = os.path.join("intermediate-files", "01-setup",\
                          "02-trim")
        extensionname = "trimmed"
    if filterrrnas == 1:
        foldername = os.path.join("intermediate-files", "01-setup",\
                          "04a-ribo")
        extensionname = "ribodetector_filt"
    if leftorright == "left":
        samplelist = [os.path.join(OUTPUTDIR, foldername, sample + "_1." + extensionname + ".fastq.gz") 
                      for sample in samplelist]
    elif leftorright == "right":
        samplelist = [os.path.join(OUTPUTDIR, foldername, sample + "_2." + extensionname + ".fastq.gz") 
                      for sample in samplelist]
    else:
        samplelist_out = []
        [samplelist_out.extend(["-1 " + os.path.join(OUTPUTDIR, foldername, sample + "_1." + extensionname + ".fastq.gz"),
                            "-2 " + os.path.join(OUTPUTDIR, foldername, sample + "_2." + extensionname + ".fastq.gz")])
                      for sample in samplelist]
        samplelist = samplelist_out
    if flag != "":
        return "--pe1" + " --pe1".join(samplelist)
    else:
        return samplelist

def get_samples_commas_spades(assemblygroup, dropspike, filterrrnas, leftorright, commas = False, retlen=False,retfirst=False,retlast=False):   
    samplelist = list(SAMPLEINFO.loc[SAMPLEINFO['AssemblyGroup'] == assemblygroup]['SampleID']) 
    foldername = os.path.join("intermediate-files", "01-setup",\
                          "03-alignment-spike")
    extensionname = "clean"
    if dropspike == 0:
        foldername = os.path.join("intermediate-files", "01-setup",\
                          "02-trim")
        extensionname = "trimmed"
    
    if filterrrnas == 1:
        foldername = os.path.join("intermediate-files", "01-setup",\
                          "04a-ribo")
        extensionname = "ribodetector_filt"
        
    if leftorright == "left":
        samplelist = [os.path.join(OUTPUTDIR, foldername, sample + "_1." + extensionname + ".fastq.gz") 
                      for sample in samplelist]
    else:
        samplelist = [os.path.join(OUTPUTDIR, foldername, sample + "_2." + extensionname + ".fastq.gz") 
                      for sample in samplelist]
    if retfirst:
        return samplelist[0]
    if retlast:
        return samplelist[-1]
    if retlen:
        return len(samplelist)
    if commas:
        trailing_num=2
        if leftorright == "left":
            trailing_num=1
        list_for_spades=["--pe"+str(number_in_row)+"-"+str(trailing_num)+" "+str(filename_in)\
                         for filename_in,number_in_row in zip(samplelist,
                                                              range(1,len(samplelist)+1))]
        return " ".join(list_for_spades)
    else:
        return samplelist
    
       
# This module needs to grab all of the list of the individual files associated with the specified
# assembly group, after the scripts/make-assembly-file.py script builds said assembly groups 
# according to user specifications.  
rule rnaspades:
    input:
        left = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, REMOVERRNA,
                                                          "left", commas = False),
        right = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, REMOVERRNA,
                                                           "right", commas = False)
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05d-rnaspades",\
                     "rna_{assembly}", "transcripts.fasta")
    params:
        continue_flag = CONTINUEFLAG,
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05d-rnaspades", "rna_{assembly}"),
        fullstring = lambda filename: get_samples_multiflag_spades(filename.assembly, DROPSPIKE, REMOVERRNA, "both"),
        numsamps = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, REMOVERRNA, 
                                                              "right", commas = False, retlen=True),
        left = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, REMOVERRNA, 
                                                           "left", commas = True, retfirst=False),
        right = lambda filename: get_samples_commas_spades(filename.assembly, DROPSPIKE, REMOVERRNA, 
                                                            "right", commas = True,
                                        retfirst=False),
        maxmem = MAXMEMORY,
        CPUs = MAXCPUSPERTASK * MAXTASKS
    log:
        err = os.path.join(OUTPUTDIR, "logs",\
                           "05-assembly", "05d-rnaspades", "{assembly}",\
                           "outputlog_{assembly}_merge.err"),
        out = os.path.join(OUTPUTDIR, "logs",\
                            "05-assembly", "05d-rnaspades", "{assembly}",\
                            "outputlog_{assembly}_merge.out") 
    conda: os.path.join("..", "..", "envs", "02-assembly-env.yaml")
    shell:
        '''
        echo {params.left}
        if [ -f {params.outdir}/params.txt ] && [ {params.continue_flag} ]; then
            spades.py --restart-from last -m 500 -o {params.outdir} 2> {log.err} 1> {log.out}
        else
            spades.py -m 100 -t 8 --rna {params.right} {params.left} -o {params.outdir} 2> {log.err} 1> {log.out}
        fi
        '''
   
rule rnaspades_cleanup:
    input:
        spadesfile = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly",\
                     "05-assembly", "05d-rnaspades",\
                     "rna_{assembly}", "transcripts.fasta")
    output:
        assembled = os.path.join(ASSEMBLEDDIR, "{assembly}_rnaspades.fasta")
    params:
        extra = ""
    shell:
        '''
        cp {input.spadesfile} {output.assembled}
        '''
