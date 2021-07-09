configfile: "config.yaml"

import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

rule transdecoder_CAG:
    input:
        fastafile = os.path.join(OUTPUTDIR, "03-merge", "07-CAG",\
                                 "{folder}", "{assembly}_merged.fasta")
    output:
        pep = os.path.join("{assembly}_CAG.fasta.transdecoder.pep"),
        gff = os.path.join("{assembly}_CAG.fasta.transdecoder.gff3"),
        cds = os.path.join("{assembly}_CAG.fasta.transdecoder.cds"),
        bed = os.path.join("{assembly}_CAG.fasta.transdecoder.bed")
    params:
        merged = "{assembly}_CAG",
        size = TRANSDECODERORFSIZE
    log:
        err = os.path.join(OUTPUTDIR, "logs", "transdecoder", "{assembly}.err"),
        out = os.path.join(OUTPUTDIR, "logs", "transdecoder", "{assembly}.log")
    conda: 
        os.path.join("..", "envs", "05-compare-env.yaml")
    shell:
        """
        unset PERL5LIB
        cp {input.merged} {params.merged}.fasta
        TransDecoder.LongOrfs -t {params.merged}.fasta -m {params.size} 2> {log.err} 1> {log.out}
        TransDecoder.Predict -t {params.merged}.fasta --no_refine_starts 2>> {log.err} 1>> {log.out}
        rm {params.merged}.fasta
        """
        
rule transdecoder_CAG_clean:
    input:
        pep = os.path.join("{assembly}_CAG.fasta.transdecoder.pep"),
        gff = os.path.join("{assembly}_CAG.fasta.transdecoder.gff3"),
        cds = os.path.join("{assembly}_CAG.fasta.transdecoder.cds"),
        bed = os.path.join("{assembly}_CAG.fasta.transdecoder.bed")
    output:
        pep = os.path.join(OUTPUTDIR, "04-compare", "08-CAG-proteins", "{assembly}.fasta.transdecoder.pep"),
        gff = os.path.join(OUTPUTDIR, "04-compare", "08-CAG-proteins", "{assembly}.fasta.transdecoder.gff3"),
        cds = os.path.join(OUTPUTDIR, "04-compare", "08-CAG-proteins", "{assembly}.fasta.transdecoder.cds"),
        bed = os.path.join(OUTPUTDIR, "04-compare", "08-CAG-proteins", "{assembly}.fasta.transdecoder.bed")
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