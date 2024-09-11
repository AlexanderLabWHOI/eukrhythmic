configfile: "config.yaml"

import io
import os
import pathlib
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

KEGG_PROT_DB = config["kegg_prot_db"]
KEGG_PATH = config["kegg"]
BUSCO_DATABASES = list(config["busco"])
PFAM = config["pfam"]

download_data=False
if not os.path.isfile(os.path.join(EGGNOG_DATA_LOC,"eggnog_proteins.dmnd")):
    download_data=True
    
rule split_files:
    input:
        os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "{filter_workflow}",\
                           "13-MAD-proteins", "MAD.fasta.transdecoder.pep")
    output:
        temp(expand(os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "{filter_workflow}", "13-MAD-proteins",
                           "MAD.fasta.transdecoder.pep.part-{splitno}"),filter_workflow="{filter_workflow}",
                           splitno=[str(curr).zfill(2) for curr in list(range(1,numbersplits))]))
    params:
        td_dir = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "{filter_workflow}",\
                           "13-MAD-proteins")
    shell:
        """
        c_wd=$(pwd)
        (cd {params.td_dir} && perl $c_wd/fasta-splitter.pl --n-parts 30 $c_wd/{input})
        """

rule emappermad_split:
    input:
        assembly_file = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "{filter_workflow}", "13-MAD-proteins", 
                           "MAD.fasta.transdecoder.pep.part-{splitno}")
    output:
        hits_file = temp(os.path.join(OUTPUTDIR, "intermediate-files",
                                 "04-compare", "{filter_workflow}",
                                 "17-MAD-emapper","split_{splitno}",
                                 "MAD_{splitno}.emapper.seed_orthologs"))
    conda: "../envs/04-compare-env.yaml"
    params:
        outdir = os.path.join(OUTPUTDIR, "intermediate-files",
                              "04-compare", "{filter_workflow}",
                              "17-MAD-emapper","split_{splitno}",),
        tmpdir = os.path.join(SCRATCHDIR,"tmp_emapper_MAD_{filter_workflow}_split_{splitno}"),
        prefix = "MAD_{splitno}",
        eggnog_mapper_data = EGGNOG_DATA_LOC,
        download_data = False
    shell:
        '''
        if [ {params.download_data} == "True" ]; then
            mkdir -p {params.eggnog_mapper_data}
            download_eggnog_data.py -y --data_dir {params.eggnog_mapper_data}
        fi
        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}
        export EGGNOG_DATA_DIR={params.eggnog_mapper_data}
        emapper.py -m diamond --override --no_annot --no_file_comments --cpu 18 -i {input.assembly_file} -o {params.prefix} --output_dir {params.outdir}
        '''

rule write_eggnog_annots:
    input:
        split_files=expand(os.path.join(OUTPUTDIR, "intermediate-files",
                                 "04-compare", "{filter_workflow}",
                                 "17-MAD-emapper","split_{splitno}",
                                 "MAD_{splitno}.emapper.seed_orthologs"),filter_workflow="{filter_workflow}",
                           splitno=[str(curr).zfill(2) for curr in list(range(1,numbersplits))])
    output:
        final_file=os.path.join(OUTPUTDIR, "intermediate-files",      
                      "04-compare", "{filter_workflow}",
                      "17-MAD-emapper","split_unite","MAD.emapper.annotations")
    params:
        outdir = os.path.join(OUTPUTDIR, "intermediate-files",
                              "04-compare","{filter_workflow}","17-MAD-emapper","split_unite"),
        tmpdir = os.path.join(SCRATCHDIR,"tmp_emapper_MAD_{filter_workflow}_split_unite"),
        prefix = "MAD",
        eggnog_mapper_data = EGGNOG_DATA_LOC,
        download_data = False
    shell:
        """
        cat {input.split_files} > {params.outdir}/MAD.emapper.seed_orthologs
        export EGGNOG_DATA_DIR={params.eggnog_mapper_data}
        emapper.py --override --annotate_hits_table {params.outdir}/MAD.emapper.seed_orthologs --no_file_comments -o {params.prefix} --output_dir {params.outdir} --cpu 18
        """
    
rule emappermad:
    input:
        assembly_file = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "{filter_workflow}",\
                           "13-MAD-proteins", "MAD.fasta.transdecoder.pep")
    output:
        hits_file = os.path.join(OUTPUTDIR, "intermediate-files",
                                 "04-compare", "{filter_workflow}", "17-MAD-emapper",
                                 "MAD.emapper.hits")
    conda: "../envs/04-compare-env.yaml"
    params:
        outdir = os.path.join(OUTPUTDIR, "intermediate-files",
                              "04-compare", "{filter_workflow}", "17-MAD-emapper"),
        tmpdir = os.path.join(SCRATCHDIR,"tmp_emapper_MAD_{filter_workflow}"),
        prefix = "MAD",
        eggnog_mapper_data = EGGNOG_DATA_LOC,
        download_data = download_data
    shell:
        '''
        if [ {params.download_data} == "True" ]; then
            mkdir -p {params.eggnog_mapper_data}
            download_eggnog_data.py -y --data_dir {params.eggnog_mapper_data}
        fi
        mkdir -p {params.outdir}
        mkdir -p {params.tmpdir}
        export EGGNOG_DATA_DIR={params.eggnog_mapper_data}
        emapper.py --override -i {input.assembly_file} --itype proteins -m diamond -o {params.prefix} --output_dir {params.outdir} --temp_dir {params.tmpdir}
        '''


