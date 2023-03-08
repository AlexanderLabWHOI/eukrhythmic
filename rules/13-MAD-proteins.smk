configfile: "config.yaml"

import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *


rule transdecoder_MAD_simple:
    input:
        fastafile = os.path.join(OUTPUTDIR, "intermediate-files",\
                           "03-merge", "{filter_workflow}", "12-MAD", "MAD.fasta")
    output:
        pep = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "{filter_workflow}",\
                           "13-MAD-proteins", "MAD.fasta.transdecoder.pep"),
        gff = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "{filter_workflow}",\
                           "13-MAD-proteins", "MAD.fasta.transdecoder.gff3"),
        cds = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare","{filter_workflow}",\
                           "13-MAD-proteins", "MAD.fasta.transdecoder.cds"),
        bed = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "{filter_workflow}",\
                           "13-MAD-proteins", "MAD.fasta.transdecoder.bed")
    params:
        merged = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "{filter_workflow}",\
                           "13-MAD-proteins","MAD"),
        filename = "MAD",
        size = TRANSDECODERORFSIZE,
        wd_path = os.path.join(OUTPUTDIR, "intermediate-files", "04-compare", "{filter_workflow}",\
                           "13-MAD-proteins")
    log:
        err = os.path.join(OUTPUTDIR, "logs", "13-MAD-proteins-{filter_workflow}", "MAD.err"),
        out = os.path.join(OUTPUTDIR, "logs", "13-MAD-proteins-{filter_workflow}", "MAD.log")
    conda: 
        os.path.join("..", "envs", "04-compare-env.yaml")
    shell:
        """
        unset PERL5LIB
        rm -rf {params.merged}*
        mkdir -p {params.wd_path}
        cp {input.fastafile} {params.merged}.fasta
        (cd {params.wd_path} && TransDecoder.LongOrfs -t {params.filename}.fasta -m {params.size}) 2> {log.err} 1> {log.out}
        (cd {params.wd_path} && TransDecoder.Predict -t {params.filename}.fasta --no_refine_starts) 2>> {log.err} 1>> {log.out}
        rm {params.merged}.fasta
        """