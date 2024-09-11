configfile: "config.yaml"

import io
import os
import pathlib
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

rule link_nt_assembly:
    input:
        os.path.join(OUTPUTDIR, "intermediate-files",
                     "03-merge", "{workflow}", "12-MAD", "MAD.fasta")
    output:
        os.path.join(OUTPUTDIR, "final-files", "00-nucleotide_assembly", 
                     "MAD.{workflow}.fasta")
    shell:
        """
        ln -s {input} {output}
        """
        
rule link_pep_assembly:
    input:
        pep = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "{workflow}",\
                           "13-MAD-proteins", "MAD.fasta.transdecoder.pep"),
        cds = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare","{workflow}",\
                           "13-MAD-proteins", "MAD.fasta.transdecoder.cds")
    output:
        pep = os.path.join(OUTPUTDIR, "final-files", "01-predicted_proteins", 
                           "MAD.{workflow}.pep"),
        cds = os.path.join(OUTPUTDIR, "final-files", "01-predicted_proteins", 
                           "MAD.{workflow}.cds")
    shell:
        """
        ln -s {input.pep} {output.pep}
        ln -s {input.cds} {output.cds}
        """
        
rule link_salmon_files:
    input:
        salmon_files = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", 
                                    "{workflow}", "14-MAD-mapping",
                                    "salmon_sample", "{assembly}_quant",
                                    "quant.sf")
    output:
        salmon_files = os.path.join(OUTPUTDIR, "final-files", "03-abundance_tables", 
                                    "salmon_{workflow}","{assembly}_quant",
                                    "quant.sf")
    shell:
        """
        ln -s {input.salmon_files} {output.salmon_files}
        """

rule link_concord_table:
    input:
        concord = os.path.join(OUTPUTDIR, "MAD.{workflow}.concordance.tsv")
    output:
        concord = os.path.join(OUTPUTDIR, "final-files", "concordance.{workflow}.tsv")
    shell:
        """
        ln -s {input.concord} {output.concord}
        """

rule link_taxonomy_and_function:
    input:
        function_file = os.path.join(OUTPUTDIR, "intermediate-files",      
                                     "04-compare", "{workflow}",
                                     "17-MAD-emapper","split_unite",
                                     "MAD.emapper.annotations"),
        taxonomy_file = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",
                                     "{workflow}",
                                     "16-MAD-taxonomy","taxonomy_estimation",
                                     "MAD.fasta.transdecoder-estimated-taxonomy.out")
    output:
        taxonomy_table = os.path.join(OUTPUTDIR, "final-files", "02-annotation_table", 
                                      "{workflow}","full_taxonomy.csv"),
        function_table = os.path.join(OUTPUTDIR, "final-files", "02-annotation_table", 
                                      "{workflow}","full_function.csv")
    shell:
        """
        ln -s {input.function_file} {output.function_table}
        ln -s {input.taxonomy_file} {output.taxonomy_table}
        """
        
rule create_funct_tax_table:
    input:
        function_file = os.path.join(OUTPUTDIR, "intermediate-files",      
                                     "04-compare", "{workflow}",
                                     "17-MAD-emapper","split_unite",
                                     "MAD.emapper.annotations"),
        taxonomy_file = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",
                                     "{workflow}",
                                     "16-MAD-taxonomy","taxonomy_estimation",
                                     "MAD.fasta.transdecoder-estimated-taxonomy.out")
    output:
        combined_table = os.path.join(OUTPUTDIR, "final-files", "02-annotation_table", 
                                      "{workflow}","TaxonomicAndFunctionalAnnotations.csv")
    run:
        function_file=pd.read_csv(input.function_file,sep="\t",header=None,comment="#",
                                  names=["query","seed_ortholog","evalue",
                                         "score","eggNOG_OGs","max_annot_lvl",
                                         "COG_category","Description","Preferred_name",
                                         "GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module",
                                         "KEGG_Reaction","KEGG_rclass","BRITE",
                                         "KEGG_TC","CAZy","BiGG_Reaction","PFAMs"])
        taxonomy_file=pd.read_csv(input.taxonomy_file,sep="\t")
        combined_file=function_file[["query","seed_ortholog","COG_category",
                                     "GOs","KEGG_ko","PFAMs"]].merge(taxonomy_file,left_on="query",
                                                                     right_on="transcript_name").\
            rename({"query":"SequenceID"},axis="columns")
        combined_file.to_csv(output.combined_table,sep="\t")
        

rule create_funct_tax_table_cag:
    input:
        function_file = expand(os.path.join(OUTPUTDIR, "intermediate-files",
                                 "04-compare", "19-CAG-emapper",
                                 "{assembly}.emapper.hits"),
                                 assembly=assemblygroups),
        taxonomy_file = expand(os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",
                                 "18-CAG-taxonomy","taxonomy_estimation",
                                 "{assembly}_CAG.fasta.transdecoder-estimated-taxonomy.out"),
                                 assembly=assemblygroups)
    output:
        combined_table = os.path.join(OUTPUTDIR, "final-files", "02-annotation_table",
            "TaxonomicAndFunctionalAnnotations_CAG.csv")
    run:
        function_files=pd.DataFrame()
        for function_file_path in input.function_file:
            function_file=pd.read_csv(function_file_path,sep="\t",header=None,comment="#",
                                      names=["query","seed_ortholog","evalue",
                                             "score","eggNOG_OGs","max_annot_lvl",
                                             "COG_category","Description","Preferred_name",
                                             "GOs","EC","KEGG_ko","KEGG_Pathway","KEGG_Module",
                                             "KEGG_Reaction","KEGG_rclass","BRITE",
                                             "KEGG_TC","CAZy","BiGG_Reaction","PFAMs"])
            function_file["AssemblyGroup"] = function_file_path.split("/")[-1].split(".emapper")[0]
            function_files=pd.concat([function_files,function_file])
        taxonomy_files=pd.DataFrame()
        for taxonomy_file_path in input.taxonomy_file:
            taxonomy_file=pd.read_csv(taxonomy_file_path,sep="\t")
            taxonomy_file["AssemblyGroup"] = taxonomy_file_path.split("/")[-1].split("_CAG")[0]
            taxonomy_files=pd.concat([taxonomy_files,taxonomy_file])
        combined_file=function_files[["query","seed_ortholog","COG_category",
                                     "GOs","KEGG_ko","PFAMs"]].merge(taxonomy_files,left_on="query",
                                                                     right_on="transcript_name").\
            rename({"query":"SequenceID"},axis="columns")
        combined_file.to_csv(output.combined_table,sep="\t")