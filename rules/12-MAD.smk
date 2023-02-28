configfile: "config.yaml"

import io
import os
from snakemake.exceptions import print_exception, WorkflowError
import sys
from Bio import SeqIO
sys.path.insert(1, '../scripts')
from importworkspace import *

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
    
## create concordance between all IDs in MAD to be used in place of cumbersome name
rule create_id_concord:
    input:
        mad = os.path.join(OUTPUTDIR, "intermediate-files",\
                           "03-merge", "12-MAD-intermed", "MAD.{filter}.fasta"),
    output:
        concordance = os.path.join(OUTPUTDIR,"MAD.{filter}.concordance.tsv")
    run:
        mad_ids = [curr.id for curr in SeqIO.parse(input.mad,"fasta")]
        pd.DataFrame({"MAD_id":mad_ids,"Concordance_Short_ID":["Sequence_"+str(curr) for curr in list(range(len(mad_ids)))]}).to_csv(output.concordance,sep=" ")
        
rule replace_mad_ids:
    input:
        mad = os.path.join(OUTPUTDIR, "intermediate-files",\
                           "03-merge", "12-MAD-intermed", "MAD.{filter}.fasta"),
        concordance = os.path.join(OUTPUTDIR,"MAD.{filter}.concordance.tsv")
    output:
        mad = os.path.join(OUTPUTDIR, "intermediate-files",\
                           "03-merge", "12-MAD", "MAD.{filter}.fasta")
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
