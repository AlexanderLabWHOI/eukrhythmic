configfile: "config.yaml"

import io
import os
from snakemake.exceptions import print_exception, WorkflowError
import sys
from Bio import SeqIO
sys.path.insert(1, '../scripts')
from importworkspace import *

def select_right_MAD(workflow):
    if workflow == "madfiltered":
        return os.path.join(OUTPUTDIR,"intermediate-files", "03-merge",
                            "merged", "20-MAD-filtered", "MAD.filtered.fasta")
    elif workflow == "fullfiltered":
        return os.path.join(OUTPUTDIR,"intermediate-files", "03-merge",
                            "filtered", "20-MAD-filtered", "MAD.filtered.fasta")
    elif workflow == "cagfiltered":
        return  os.path.join(OUTPUTDIR, "intermediate-files",\
                             "03-merge", "12-MAD-intermed", "MAD.filtered.fasta")
    elif workflow == "nofilter":
        return  os.path.join(OUTPUTDIR, "intermediate-files",\
                             "03-merge", "12-MAD-intermed", "MAD.merged.fasta")

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

def get_samples(assemblygroup):
    samplelist = list(ASSEMBLYFILE.loc[ASSEMBLYFILE['AssemblyGroup'] == assemblygroup]['SampleID']) 
    return samplelist

rule mad_mmseqs:
    input: 
        infiles = os.path.join(OUTPUTDIR, "intermediate-files", "02-assembly", "11-SWAM", "{filter}.fasta")
    output:
        outfasta = os.path.join(OUTPUTDIR, "intermediate-files",\
                                "03-merge", "12-MAD-intermed", "MAD.{filter}.fasta"),
        outtsv = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                              "12-MAD-intermed", "MAD.{filter}.tsv")
    params:
        threads = 10,
        maxmemory = 30000, # -G o indicates local sequence identity.
        identityparam = 1.00,
        mincoverageshorter = MINCOVERAGECLUST2,
        mincoveragelong = 0.005,
        name_db = "MAD_mmseqs",
        name_intermed = "MAD_mmseqs_2",
        name_subdb = "MAD_mmseqs_3"
    log:
        err = os.path.join(OUTPUTDIR, "logs", "12-MAD", "MAD.{filter}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "12-MAD", "MAD.{filter}.log")
    conda:
        os.path.join("..", "envs", "03-merge-env.yaml")
    shell:
        '''
        mmseqs createdb {input.infiles} {params.name_db} 
        mmseqs linclust {params.name_db} {params.name_intermed} tmp --min-seq-id {params.identityparam} --cov-mode 1 -c {params.mincoverageshorter} --split-memory-limit 120G --remove-tmp-files 2> {log.err} 1> {log.out}
        mmseqs createsubdb {params.name_intermed} {params.name_db} {params.name_subdb}
        mmseqs convert2fasta {params.name_subdb} {output.outfasta}
        mmseqs createtsv {params.name_db} {params.name_db} {params.name_intermed} {output.outtsv}
        rm -f {params.name_db}*
        rm -f {params.name_intermed}*
        rm -f {params.name_subdb}*
        '''

## select the file that we want to be used as our final MAD. this could be either an unfiltered file,
## a file filtered by CAG only, MAD only, or both CAG+MAD
rule select_proper_file:
    input:
        selected_file = lambda filename: select_right_MAD(filename.filter_workflow)
    output:
        new_mad_loc = os.path.join(OUTPUTDIR, "intermediate-files",
                           "03-merge", "{filter_workflow}", "12-MAD", "MAD.fasta")
    shell:
        """
        mv {input.selected_file} {output.new_mad_loc}
        """
        
## create concordance between all IDs in MAD to be used in place of cumbersome name
rule create_id_concord:
    input:
        mad =  os.path.join(OUTPUTDIR, "intermediate-files",\
                            "03-merge", "12-MAD-intermed", "MAD.{filter}.fasta")
    output:
        concordance = os.path.join(OUTPUTDIR,"MAD.{filter_workflow}.concordance.tsv")
    run:
        mad_ids = [curr.id for curr in SeqIO.parse(input.mad,"fasta")]
        pd.DataFrame({"MAD_id":mad_ids,"Concordance_Short_ID":["Sequence_"+str(curr) for curr in list(range(len(mad_ids)))]}).to_csv(output.concordance,sep=" ")
        
rule replace_mad_ids:
    input:
        mad = os.path.join(OUTPUTDIR, "intermediate-files",\
                           "03-merge", "12-MAD-intermed", "MAD.{filter}.fasta"),
        concordance = os.path.join(OUTPUTDIR,"MAD.{filter}.concordance.tsv")
    output:
        mad = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "12-MAD-intermed",\
                           "rename","MAD.{filter}.fasta")
    shell:
        """
        ## modified from https://unix.stackexchange.com/questions/652523/replacing-the-seq-ids-of-fasta-file-based-on-the-new-ids-from-a-list
        awk -F',' '
          NR==FNR{{ a[$1]=$2; next }}
          /^>/{{ 
            id=a[substr($0, 2)]
            if (id!=""){{ print ">" id; next }}
          }}
          1
        ' {input.concordance} {input.mad} > {output.mad}
        """

rule convert_mad_no_space_temp:
    input:
        fastafile = os.path.join(OUTPUTDIR, "intermediate-files",
                                 "03-merge", "{filter_workflow}", "12-MAD-intermed", "MAD.fasta")
        #os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "12-MAD",\
        #                         "MAD.{filter}.fasta")
    output:
        fastafile = temp(os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "12-MAD",\
                                 "MAD.{filter}.nospace.fasta"))
    shell:
        '''
        sed 's, ,_,g' {input.fastafile} > {output.fastafile}
        '''

rule convert_mad_no_space_temp_intermed:
    input:
        fastafile = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "12-MAD-intermed",\
                                 "rename","MAD.{filter}.fasta")
    output:
        fastafile = temp(os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "12-MAD-intermed",\
                                 "rename","MAD.{filter}.nospace.fasta"))
    shell:
        '''
        sed 's, ,_,g' {input.fastafile} > {output.fastafile}
        '''
        
rule salmon_MAD_index_temp:
    input: 
        fastafile = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "12-MAD-intermed",\
                                 "rename","MAD.{filter}.nospace.fasta")
    output:
        temp(os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping",\
                     "salmon", "MAD_{filter}_index", "versionInfo.json"))
    params:
        libtype = "A",
        indexname = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping",\
                     "salmon", "MAD_{filter}_index"),
        kval = 31
    log:
        err = os.path.join(OUTPUTDIR, "logs", "14-MAD-mapping", "MAD_salmon_ind_{filter}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "14-MAD-mapping", "MAD_salmon_ind_{filter}.out")
    conda: os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        """
        salmon index -t {input.fastafile} -i {params.indexname} -k {params.kval} 2> {log.err} 1> {log.out}
        """
        
rule salmon_MAD_temp:
    input: 
        indexname = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping",\
                     "salmon", "MAD_{filter}_index", "versionInfo.json"),
        left = lambda filename: salmon_get_sample(filename.assembly, "left", list_format = True),
        right = lambda filename: salmon_get_sample(filename.assembly, "right", list_format = True)
    output:
        temp(os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping",\
                     "salmon_sample_{filter}", "{assembly}_quant", "quant.sf"))
    params:
        libtype = "A",
        indexname = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping",\
                     "salmon_sample_{filter}", "MAD_{filter}_index"),
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping",\
                     "salmon_sample_{filter}", "{assembly}_quant"),
        kval = 31
    log:
        err = os.path.join(OUTPUTDIR, "logs", "14-MAD-mapping", "{assembly}_{filter}_salmon.err"),
        out = os.path.join(OUTPUTDIR, "logs", "14-MAD-mapping", "{assembly}_{filter}_salmon.out")
    conda: os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        """
        salmon quant -i {params.indexname} -l {params.libtype} -1 {input.left} -2 {input.right} -p 20 --validateMappings -o {params.outdir} 2> {log.err} 1> {log.out}
        """
        
rule mad_filter_by_salmon:
    input:
        salmon_files = expand(os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", 
                                           "14-MAD-mapping",\
                     "salmon_sample_{filter}", "{assembly}_quant", "quant.sf"),assembly=filenames,filter="{filter}"),
        MAD_file = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",
                                 "12-MAD-intermed",\
                                 "rename","MAD.{filter}.nospace.fasta")
    output:
        seq_extract = temp(os.path.join(OUTPUTDIR,"intermediate-files",
                                        "03-merge", "{filter}", "20-MAD-filtered", "MAD.seqs.txt"))
    run:
        all_quant=pd.DataFrame()
        for quant_file in input.salmon_files:
            if (not "quant" in quant_file)|("merged" in quant_file):
                continue
            current_quant=pd.read_csv(os.path.join(quant_file),sep="\t")
            current_quant=current_quant.loc[current_quant.NumReads>0,["Name","NumReads"]]
            current_quant["Sample"] = quant_file.split("_quant")[0].split("/")[0]
            all_quant=pd.concat([all_quant,current_quant])
            
        summed_quant=all_quant.groupby(["Name"]).NumReads.sum().reset_index()
        zerod_names=summed_quant.loc[summed_quant.NumReads==0,"Name"]
        to_keep_seq=list(set(summed_quant.Name)-set(zerod_names))
        with open(params.seq_extract,"w") as f:
            f.write("\n".join(to_keep_seq))

            
rule seqtk_by_salmon:
    input:
        seq_extract = os.path.join(OUTPUTDIR,"intermediate-files", "03-merge", "{filter}", 
                                   "20-MAD-filtered", "MAD.seqs.txt"),
        salmon_files = expand(os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", 
                                           "14-MAD-mapping",\
                     "salmon_sample_{filter}", "{assembly}_quant", "quant.sf"),assembly=filenames,filter="{filter}"),
        MAD_file = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "12-MAD-intermed",\
                                 "rename","MAD.{filter}.nospace.fasta")
    output:
        os.path.join(OUTPUTDIR,"intermediate-files", "03-merge", "{filter}", "20-MAD-filtered", "MAD.filtered.fasta")
    conda: "../envs/03-merge-env.yaml"
    shell:
        """
        seqtk subseq {input.MAD_file} {input.seq_extract} > {output}
        """
        
