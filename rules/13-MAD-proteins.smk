configfile: "config.yaml"

import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

rule transdecoder_MAD:
    input:
        fastafile = os.path.join(OUTPUTDIR, "intermediate-files", "03-merge", "12-MAD", "MAD.fasta")
    output:
        pep = os.path.join("MAD.fasta.transdecoder.pep"),
        gff = os.path.join("MAD.fasta.transdecoder.gff3"),
        cds = os.path.join("MAD.fasta.transdecoder.cds"),
        bed = os.path.join("MAD.fasta.transdecoder.bed")
    params:
        merged = "TD_MAD",
        size = TRANSDECODERORFSIZE
    log:
        err = os.path.join(OUTPUTDIR, "logs", "13-MAD-proteins", "MAD.err"),
        out = os.path.join(OUTPUTDIR, "logs", "13-MAD-proteins", "MAD.log")
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
        
rule transdecoder_MAD_clean:
    input:
        pep = os.path.join("MAD.fasta.transdecoder.pep"),
        gff = os.path.join("MAD.fasta.transdecoder.gff3"),
        cds = os.path.join("MAD.fasta.transdecoder.cds"),
        bed = os.path.join("MAD.fasta.transdecoder.bed")
    output:
        pep = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "13-MAD-proteins", "MAD.fasta.transdecoder.pep"),
        gff = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "13-MAD-proteins", "MAD.fasta.transdecoder.gff3"),
        cds = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "13-MAD-proteins", "MAD.fasta.transdecoder.cds"),
        bed = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare",\
                           "13-MAD-proteins", "MAD.fasta.transdecoder.bed")
    params:
        merged = "TD_MAD",
        size = TRANSDECODERORFSIZE
    shell:
        """
        mv {input.pep} {output.pep}
        mv {input.cds} {output.cds}
        mv {input.gff} {output.gff}
        mv {input.bed} {output.bed}
        rm -rf {params.merged}.fasta.transdecoder_dir*
        rm -rf pipeliner.*.cmds
        """