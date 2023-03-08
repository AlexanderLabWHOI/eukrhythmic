configfile: "config.yaml"

import io
import os
import sys
from snakemake.exceptions import print_exception, WorkflowError
import sys
import subprocess
sys.path.insert(1, '../scripts')
from importworkspace import *

download_data=False
if not os.path.isfile(os.path.join(EGGNOG_DATA_LOC,"eggnog_proteins.dmnd")):
    download_data=True
sep_list=["0" + str(curr) if curr < 10 else str(curr) for curr in list(range(50))]
    
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
        
rule cag_filter_by_salmon:
    input:
        salmon_files = expand(os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "09-CAG-mapping",\
                     "salmon", "{assembly}_quant", "quant.sf"),assembly=assemblygroups),
        CAG_files = expand(os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "07-CAG",\
                                 "{assembly}_merged.fasta"),assembly=assemblygroups)
    output:
        seq_extract = temp(expand(os.path.join(OUTPUTDIR,"intermediate-files", "03-merge",
                                               "19-CAG-filtered", "{assembly}.seqs.txt"),
                                  assembly=assemblygroups))
    run:
        all_quant=pd.DataFrame()
        for quant_file,cag_file,out_file in zip(input.salmon_files,input.CAG_files,output.seq_extract):
            if (not "quant" in quant_file)|("merged" in quant_file):
                continue
            current_quant=pd.read_csv(os.path.join(quant_file),sep="\t")
            current_quant=current_quant.loc[current_quant.NumReads>0,["Name","NumReads"]]
            current_quant["Sample"] = quant_file.split("_quant")[0].split("/")[0]
            all_quant=pd.concat([all_quant,current_quant])
            
            summed_quant=all_quant.groupby(["Name"]).NumReads.sum().reset_index()
            zerod_names=summed_quant.loc[summed_quant.NumReads==0,"Name"]
            to_keep_seq=list(set(summed_quant.Name)-set(zerod_names))
            with open(out_file,"w") as f:
                f.write("\n".join(to_keep_seq))
                
            
rule seqtk_by_salmon_cag:
    input:
        seq_extract = temp(os.path.join(OUTPUTDIR,"intermediate-files", "03-merge",
                                               "19-CAG-filtered", "{assembly}.seqs.txt")),
        salmon_files = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "09-CAG-mapping",\
                     "salmon", "{assembly}_quant", "quant.sf"),
        CAG_file = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "07-CAG",\
                                 "{assembly}_merged.fasta")
    output:
        os.path.join(OUTPUTDIR,"intermediate-files", "03-merge", "19-CAG-filtered",
                            "{assembly}.filtered.fasta")
    conda: "../envs/03-merge-env.yaml"
    shell:
        """
        seqtk subseq {input.CAG_file} {input.seq_extract} > {output}
        """
        
rule salmon_against_filtered_index_cag:
    input:
        filtered_cag_file = os.path.join(OUTPUTDIR,"intermediate-files", "03-merge", "19-CAG-filtered",
                            "{assembly}.filtered.fasta")
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "19-CAG-filtered",\
                     "salmon", "{assembly}_index", "versionInfo.json")
    params:
        libtype = "A",
        indexname = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "19-CAG-filtered",\
                     "salmon", "{assembly}_index"),
        kval = 31
    log:
        err = os.path.join(OUTPUTDIR, "logs", "19-CAG-filtered", "{assembly}_salmon_ind.err"),
        out = os.path.join(OUTPUTDIR, "logs", "19-CAG-filtered", "{assembly}_salmon_ind.out")
    conda: os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        """
        salmon index -t {input.filtered_mad_file} -i {params.indexname} -k {params.kval} 2> {log.err} 1> {log.out}
        """
        
rule salmon_against_filtered_cag:
    input: 
        indexname = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "19-CAG-filtered",\
                     "salmon", "{assembly}_index", "versionInfo.json"),
        left = lambda filename: salmon_get_sample(filename.sample, "left", list_format = True),
        right = lambda filename: salmon_get_sample(filename.sample, "right", list_format = True)
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "19-CAG-filtered",\
                     "salmon_sample", "{assembly}_quant", "sample_{sample}" "quant.sf")
    params:
        libtype = "A",
        indexname = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "19-CAG-filtered",\
                     "salmon", "{assembly}_index"),
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "19-CAG-filtered",\
                     "salmon_{assembly}_sample", "{sample}_quant"),
        kval = 31
    log:
        err = os.path.join(OUTPUTDIR, "logs", "19-CAG-filtered","{assembly}_group","{sample}_salmon.err"),
        out = os.path.join(OUTPUTDIR, "logs", "19-CAG-filtered","{assembly}_group","{sample}_salmon.out")
    conda: os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        """
        salmon quant -i {params.indexname} -l {params.libtype} -1 {input.left} -2 {input.right} -p 20 --validateMappings -o {params.outdir} 2> {log.err} 1> {log.out}
        """