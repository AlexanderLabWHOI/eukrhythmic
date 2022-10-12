configfile: "config.yaml"

import io
import os
import sys
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *
    
def salmon_get_samples(assembly,left_or_right,list_format):
    foldername = os.path.join("intermediate-files", "01-setup",\
                          "03-alignment-spike")
    extensionname = "clean"
    if DROPSPIKE == 0:
        foldername = os.path.join("intermediate-files", "01-setup",\
                          "02-trim")
        extensionname = "trimmed"
    samplelist = list(SAMPLEINFO.loc[SAMPLEINFO['AssemblyGroup'] == assembly]['SampleID']) 
    if assembly == "merged":
        samplelist = list(SAMPLEINFO['SampleID'])

    if left_or_right == "left":
        filenames = [os.path.join(OUTPUTDIR, foldername, sample + "_1." + extensionname + ".fastq.gz") 
                      for sample in samplelist]
    else:
        filenames = [os.path.join(OUTPUTDIR, foldername, sample + "_2." + extensionname + ".fastq.gz") 
                      for sample in samplelist]
    if list_format:
        return filenames
    else:
        return " ".join(filenames)
        

rule salmon_CAG:
    input: 
        fastafile = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "07-CAG",\
                                 "{assembly}_merged.fasta"),
        left = lambda filename: salmon_get_samples(filename.assembly, "left", list_format = True),
        right = lambda filename: salmon_get_samples(filename.assembly, "right", list_format = True)
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "09-CAG-mapping",\
                     "salmon", "{assembly}_quant", "quant.sf")
    params:
        libtype = "A",
        indexname = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "09-CAG-mapping",\
                     "salmon", "{assembly}_index"),
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "09-CAG-mapping",\
                     "salmon", "{assembly}_quant"),
        kval = 31
    log:
        err = os.path.join(OUTPUTDIR, "logs", "09-CAG-mapping", "{assembly}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "09-CAG-mapping", "{assembly}.log")
    conda: os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        """
        salmon index -t {input.fastafile} -i {params.indexname} -k {params.kval} 2> {log.err} 1> {log.out}
        salmon quant -i {params.indexname} -l {params.libtype} -1 {input.left} -2 {input.right} --validateMappings -o {params.outdir} 2>> {log.err} 1>> {log.out}
        """