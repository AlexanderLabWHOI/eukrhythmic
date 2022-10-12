configfile: "config.yaml"

import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

rule transdecoder_CAG_simple:
    input:
        fastafile = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "07-CAG",\
                                 "{assembly}_merged.fasta")
    output:
        pep = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "08-CAG-proteins", "{assembly}_CAG.fasta.transdecoder.pep"),
        gff = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "08-CAG-proteins", "{assembly}_CAG.fasta.transdecoder.gff3"),
        cds = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "08-CAG-proteins", "{assembly}_CAG.fasta.transdecoder.cds"),
        bed = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "08-CAG-proteins", "{assembly}_CAG.fasta.transdecoder.bed")
    params:
        merged = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "08-CAG-proteins","{assembly}_CAG"),
        filename = "{assembly}_CAG",
        wd_path = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "08-CAG-proteins"),
        size = TRANSDECODERORFSIZE
    log:
        err = os.path.join(OUTPUTDIR, "logs", "08-CAG-proteins", "{assembly}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "08-CAG-proteins", "{assembly}.log")
    conda: 
        os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        """
        unset PERL5LIB
        mkdir -p {params.wd_path}
        cp {input.fastafile} {params.merged}.fasta
        (cd {params.wd_path} && TransDecoder.LongOrfs -t {params.filename}.fasta -m {params.size}) 2> {log.err} 1> {log.out}
        (cd {params.wd_path} && TransDecoder.Predict --no_refine_starts -t {params.filename}.fasta --no_refine_starts) 2>> {log.err} 1>> {log.out}
        rm {params.merged}.fasta
        sleep 10
        """

rule cp_td:
    input:
        pep = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "08-CAG-proteins", "{assembly}_CAG.fasta.transdecoder.pep")
    output:
        pep = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "08-CAG-proteins", "eukulele_in",
                           "{assembly}_CAG.pep")
    shell:
        """
        cp {input.pep} {output.pep}
        """


rule transdecoder_CAG:
    input:
        fastafile = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "07-CAG",\
                                 "{assembly}_merged.fasta")
    output:
        pep = os.path.join("old.{assembly}_CAG.fasta.transdecoder.pep"),
        gff = os.path.join("old.{assembly}_CAG.fasta.transdecoder.gff3"),
        cds = os.path.join("old.{assembly}_CAG.fasta.transdecoder.cds"),
        bed = os.path.join("old.{assembly}_CAG.fasta.transdecoder.bed")
    params:
        merged = "{assembly}_CAG",
        size = TRANSDECODERORFSIZE
    log:
        err = os.path.join(OUTPUTDIR, "logs", "08-CAG-proteins", "{assembly}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "08-CAG-proteins", "{assembly}.log")
    conda: 
        os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        """
        unset PERL5LIB
        cp {input.fastafile} {params.merged}.fasta
        TransDecoder.LongOrfs -t {params.merged}.fasta -m {params.size} 2> {log.err} 1> {log.out}
        TransDecoder.Predict -t {params.merged}.fasta --no_refine_starts 2>> {log.err} 1>> {log.out}
        rm {params.merged}.fasta
        """
        
rule transdecoder_CAG_clean:
    input:
        pep = os.path.join("old.{assembly}_CAG.fasta.transdecoder.pep"),
        gff = os.path.join("old.{assembly}_CAG.fasta.transdecoder.gff3"),
        cds = os.path.join("old.{assembly}_CAG.fasta.transdecoder.cds"),
        bed = os.path.join("old.{assembly}_CAG.fasta.transdecoder.bed")
    output:
        pep = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "08-CAG-proteins", "{assembly}.old.fasta.transdecoder.pep"),
        gff = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "08-CAG-proteins", "{assembly}.old.fasta.transdecoder.gff3"),
        cds = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "08-CAG-proteins", "{assembly}.old.fasta.transdecoder.cds"),
        bed = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "08-CAG-proteins", "{assembly}.old.fasta.transdecoder.bed")
    params:
        merged = "{assembly}_CAG",
        size = TRANSDECODERORFSIZE
    shell:
        """
        mv {input.pep} {output.pep}
        mv {input.cds} {output.cds}
        mv {input.gff} {output.gff}
        mv {input.bed} {output.bed}
        rm -rf {params.merged}.fasta.transdecoder_dir*
        # rm -rf pipeliner.*.cmds
        """
