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
    
rule mad_filter_by_salmon:
    input:
        salmon_files = expand(os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping-merged",\
                     "salmon_sample", "{assembly}_quant", "quant.sf"),assembly=assemblygroups),
        MAD_file = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "12-MAD",\
                                 "MAD.merged.nospace.fasta")
    output:
        seq_extract = temp(os.path.join(OUTPUTDIR,"intermediate-files", "03-merge", "20-MAD-filtered", "MAD.seqs.txt"))
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
        seq_extract = os.path.join(OUTPUTDIR,"intermediate-files", "03-merge", "20-MAD-filtered", "MAD.seqs.txt"),
        salmon_files = expand(os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "14-MAD-mapping-merged",\
                     "salmon_sample", "{assembly}_quant", "quant.sf"),assembly=assemblygroups),
        MAD_file = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "12-MAD",\
                                 "MAD.merged.nospace.fasta")
    output:
        os.path.join(OUTPUTDIR,"intermediate-files", "03-merge", "20-MAD-filtered", "MAD.filtered.fasta")
    conda: "../envs/03-merge-env.yaml"
    shell:
        """
        seqtk subseq {input.MAD_file} {input.seq_extract} > {output}
        """
        
rule salmon_against_filtered_index:
    input:
        filtered_mad_file = os.path.join(OUTPUTDIR,"intermediate-files", "03-merge", "20-MAD-filtered", "MAD.filtered.fasta")
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "20-MAD-filtered",\
                     "salmon", "MAD_index", "versionInfo.json")
    params:
        libtype = "A",
        indexname = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "20-MAD-filtered",\
                     "salmon", "MAD_index"),
        kval = 31
    log:
        err = os.path.join(OUTPUTDIR, "logs", "20-MAD-filtered", "MAD_salmon_ind.err"),
        out = os.path.join(OUTPUTDIR, "logs", "20-MAD-filtered", "MAD_salmon_ind.out")
    conda: os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        """
        salmon index -t {input.filtered_mad_file} -i {params.indexname} -k {params.kval} 2> {log.err} 1> {log.out}
        """
        
rule salmon_against_filtered:
    input: 
        indexname = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "20-MAD-filtered",\
                     "salmon", "MAD_index", "versionInfo.json"),
        left = lambda filename: salmon_get_sample(filename.sample, "left", list_format = True),
        right = lambda filename: salmon_get_sample(filename.sample, "right", list_format = True)
    output:
        os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "20-MAD-filtered",\
                     "salmon_sample", "{sample}_quant", "quant.sf")
    params:
        libtype = "A",
        indexname = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "20-MAD-filtered",\
                     "salmon", "MAD_index"),
        outdir = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "20-MAD-filtered",\
                     "salmon_sample", "{sample}_quant"),
        kval = 31
    log:
        err = os.path.join(OUTPUTDIR, "logs", "20-MAD-filtered", "{sample}_salmon.err"),
        out = os.path.join(OUTPUTDIR, "logs", "20-MAD-filtered", "{sample}_salmon.out")
    conda: os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        """
        salmon quant -i {params.indexname} -l {params.libtype} -1 {input.left} -2 {input.right} -p 20 --validateMappings -o {params.outdir} 2> {log.err} 1> {log.out}
        """
       
rule separate_mad:
    input:
        fastafile = os.path.join(OUTPUTDIR,"intermediate-files", "03-merge", 
                                        "20-MAD-filtered", "MAD.filtered.fasta")
    output:
        fastafiles = expand(os.path.join(OUTPUTDIR,"intermediate-files", "03-merge", 
                                        "20-MAD-filtered", "MAD.filtered.{sep}.fasta"),sep=sep_list)
    params:
        num_splits = len(sep_list),
        location = os.path.join(OUTPUTDIR,"intermediate-files", "03-merge", 
                                        "20-MAD-filtered")
    shell:
        """
        ((cd {params.location}) && (pyfasta split -n {params.num_splits} {input.fastafile}))
        """

rule filtered_transdecoder:
    input:
        fastafile = os.path.join(OUTPUTDIR,"intermediate-files", "03-merge", 
                                        "20-MAD-filtered", "MAD.filtered.{sep}.fasta")
    output:
        pep = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                            "20-MAD-filtered", "transdecoder", "MAD.{sep}.fasta.transdecoder.pep"),
        gff = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                           "20-MAD-filtered", "transdecoder", "MAD.{sep}.fasta.transdecoder.gff3"),
        cds = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                           "20-MAD-filtered", "transdecoder", "MAD.{sep}.fasta.transdecoder.cds"),
        bed = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                           "20-MAD-filtered", "transdecoder", "MAD.{sep}.fasta.transdecoder.bed")
    params:
        merged = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                           "20-MAD-filtered", "transdecoder","sep_{sep}_sep","MAD.{sep}"),
        filename = "MAD.{sep}",
        size = TRANSDECODERORFSIZE,
        wd_path = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",
                    "20-MAD-filtered", "transdecoder","sep_{sep}_sep"),
        wd_path_top = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",
                    "20-MAD-filtered", "transdecoder")
    log:
        err = os.path.join(OUTPUTDIR, "logs", "20-MAD-filtered", "transdecoder", "MAD.{sep}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "20-MAD-filtered", "transdecoder", "MAD.{sep}.log")
    conda: 
        os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        """
        unset PERL5LIB
        mkdir -p {params.wd_path}
        cp {input.fastafile} {params.merged}.fasta
        (cd {params.wd_path} && TransDecoder.LongOrfs -t {params.filename}.fasta -m {params.size}) 2> {log.err} 1> {log.out}
        (cd {params.wd_path} && TransDecoder.Predict -t {params.filename}.fasta --no_refine_starts) 2>> {log.err} 1>> {log.out}
        rm {params.merged}.fasta
        mv {params.wd_path}/*.pep {params.wd_path_top}
        mv {params.wd_path}/*.cds {params.wd_path_top}
        mv {params.wd_path}/*.bed {params.wd_path_top}
        mv {params.wd_path}/*.gff3 {params.wd_path_top}
        rm -rf {params.wd_path}
        """
 
rule combine_transdecoder:
    input:
        pep=expand(os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                            "20-MAD-filtered", "transdecoder", "MAD.{sep}.fasta.transdecoder.pep"),sep=sep_list),
        cds=expand(os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                           "20-MAD-filtered", "transdecoder", "MAD.{sep}.fasta.transdecoder.cds"),sep=sep_list)
    output:
        pep=os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                            "20-MAD-filtered", "transdecoder", "MAD.fasta.transdecoder.pep"),
        cds=os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                           "20-MAD-filtered", "transdecoder", "MAD.fasta.transdecoder.cds")
    shell:
        """
        cat {input.pep} > {output.pep}
        cat {input.cds} > {output.cds}
        """
        
rule eukulele_run_filtered:
    input:
        pep = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                            "20-MAD-filtered", "transdecoder", "MAD.fasta.transdecoder.pep"),
        salmon = (lambda filename: expand(os.path.join(OUTPUTDIR, "intermediate-files",
                            "03-merge","20-MAD-filtered","salmon_sample","{assembly}_quant",
                            "quant.sf"), 
               assembly = filenames))
    output:
        tax_est = os.path.join(OUTPUTDIR, "intermediate-files","03-merge",
                    "20-MAD-filtered","EUKulele","taxonomy_estimation", "MAD.fasta.transdecoder-estimated-taxonomy.out")
    params:
        sampledir = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                           "20-MAD-filtered"),
        salmondir = os.path.join(OUTPUTDIR, "intermediate-files","03-merge",
                                 "20-MAD-filtered","salmon"),
        outputdir = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                           "20-MAD-filtered","EUKulele")
    #conda:
    #    os.path.join("..","envs","eukulele-env.yaml")
    shell:
        """
        EUKulele --mets_or_mags mets --sample_dir {params.sampledir} --p_ext ".pep" --n_ext ".fasta" --reference_dir /vortexfs1/omics/alexander/akrinos/remodeling/EUKulele/databases/marmmetsp_better_diatom_taxonomy --output_dir {params.outputdir}
        #--use_salmon_counts --salmon_dir {params.salmondir}
        """

rule eukulele_run_filtered_contigs:
    input:
        pep = os.path.join(OUTPUTDIR,"intermediate-files", "03-merge", 
                                        "20-MAD-filtered", "MAD.filtered.fasta"),
        salmon = (lambda filename: expand(os.path.join(OUTPUTDIR, "intermediate-files",
                            "03-merge","20-MAD-filtered","salmon_sample","{assembly}_quant",
                            "quant.sf"), 
               assembly = filenames))
    output:
        tax_est = os.path.join(OUTPUTDIR, "intermediate-files","03-merge",
                    "20-MAD-filtered","EUKulele","taxonomy_estimation", "MAD-estimated-taxonomy.out")
    params:
        sampledir = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                           "20-MAD-filtered","EUKulele"),
        salmondir = os.path.join(OUTPUTDIR, "intermediate-files","03-merge",
                                 "20-MAD-filtered","salmon")
    #conda:
    #    os.path.join("..","envs","eukulele-env.yaml")
    shell:
        """
        EUKulele --mets_or_mags mets --sample_dir {params.sampledir} --p_ext ".pep" --n_ext ".cds" --reference_dir /vortexfs1/omics/alexander/akrinos/remodeling/EUKulele/databases/marmmetsp_better_diatom_taxonomy
        #--use_salmon_counts --salmon_dir {params.salmondir}
        """
        
rule emappermad_filt:
    input:
        assembly_file = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge",\
                            "20-MAD-filtered", "transdecoder", "MAD.fasta.transdecoder.pep")
    output:
        hits_file = os.path.join(OUTPUTDIR, "intermediate-files",
                                 "03-merge", "20-MAD-filtered", "emapper",
                                 "MAD.emapper.hits")
    conda: "../envs/04-compare-env.yaml"
    params:
        outdir = os.path.join(OUTPUTDIR, "intermediate-files",
                              "03-merge", "20-MAD-filtered","emapper"),
        tmpdir = os.path.join(SCRATCHDIR,"tmp_emapper_MAD_filt"),
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
