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

## for samples, not assembly group
def salmon_get_sample(sample_in,left_or_right,list_format):
    foldername = os.path.join("intermediate-files", "01-setup",\
                          "03-alignment-spike")
    extensionname = "clean"
    if DROPSPIKE == 0:
        foldername = os.path.join("intermediate-files", "01-setup",\
                          "02-trim")
        extensionname = "trimmed"
    samplelist = sample_in

    if left_or_right == "left":
        filenames = os.path.join(OUTPUTDIR, foldername, sample_in + "_1." + extensionname + ".fastq.gz") 
    else:
        filenames = os.path.join(OUTPUTDIR, foldername, sample_in + "_2." + extensionname + ".fastq.gz")
    if list_format:
        return filenames
    else:
        return " ".join(filenames)
        
rule convert_mad_no_space:
    input:
        fastafile = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "12-MAD",\
                                 "MAD.{filter_keyword}.fasta")
    output:
        fastafile = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "12-MAD",\
                                 "MAD.{filter_keyword}.nospace.fasta")
    shell:
        '''
        sed 's, ,_,g' {input.fastafile} > {output.fastafile}
        '''
        

rule salmon_MAD_index:
    input: 
        fastafile = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "12-MAD",\
                                 "MAD.{filter_keyword}.nospace.fasta")
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping-{filter_keyword}",\
                     "salmon", "MAD.{filter_keyword}.index", "versionInfo.json")
    params:
        libtype = "A",
        indexname = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping-{filter_keyword}",\
                     "salmon", "MAD.{filter_keyword}.index"),
        kval = 31
    log:
        err = os.path.join(OUTPUTDIR, "logs", "14-MAD-mapping", "MAD_{filter_keyword}_salmon_ind.err"),
        out = os.path.join(OUTPUTDIR, "logs", "14-MAD-mapping", "MAD_{filter_keyword}_salmon_ind.out")
    conda: os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        """
        salmon index -t {input.fastafile} -i {params.indexname} -k {params.kval} 2> {log.err} 1> {log.out}
        """
        
rule salmon_MAD:
    input: 
        indexname = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping-{filter_keyword}",\
                     "salmon", "MAD.{filter_keyword}.index", "versionInfo.json"),
        left = lambda filename: salmon_get_sample(filename.sample, "left", list_format = True),
        right = lambda filename: salmon_get_sample(filename.sample, "right", list_format = True),
        gff = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "13-MAD-proteins-{filter_keyword}", "MAD.fasta.transdecoder.gff3")
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping-{filter_keyword}",\
                     "salmon_sample", "{sample}_quant", "quant.sf"),
        os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping-{filter_keyword}",\
                     "salmon_sample", "{sample}_quant", "quant.genes.sf")
    params:
        libtype = "A",
        indexname = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping-{filter_keyword}",\
                     "salmon", "MAD.{filter_keyword}.index"),
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping-{filter_keyword}",\
                     "salmon_sample", "{sample}_quant"),
        kval = 31
    log:
        err = os.path.join(OUTPUTDIR, "logs", "14-MAD-mapping", "{sample}_salmon.{filter_keyword}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "14-MAD-mapping", "{sample}_salmon.{filter_keyword}.out")
    conda: os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        """
        salmon quant -g {input.gff} -i {params.indexname} -l {params.libtype} -1 {input.left} -2 {input.right} -p 20 --validateMappings -o {params.outdir} 2> {log.err} 1> {log.out}
        """