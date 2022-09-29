configfile: "config.yaml"
  
import io
import os
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *


rule eukulele_pleuro:
    input:
        pep = (lambda filename: expand(os.path.join(OUTPUTDIR, "intermediate-files",
                         "04-compare",\
                           "08-CAG-proteins","eukulele_in",
               "{assembly}_CAG.pep"),
               assembly = assemblygroups)),
        salmon = (lambda filename: expand(os.path.join(OUTPUTDIR, "intermediate-files",\
                 "04-compare",
             "09-CAG-mapping",\
                     "salmon", "{assembly}_quant", "quant.sf"),
               assembly = assemblygroups))
    output:
        eukulele_directory = directory(os.path.join(OUTPUTDIR, "CAG_eukulele_pleuro")),
        eukulele_done = os.path.join(OUTPUTDIR, "CAG_eukulele_pleuro", "EUKulele_done.txt")
        #tax_est = os.path.join(OUTPUTDIR, "eukulele_CAG", "taxonomy_estimation", "{CAG}_)
    params:
        sampledir = os.path.join(OUTPUTDIR, "intermediate-files",
                                             "04-compare",\
                           "08-CAG-proteins","eukulele_in"),
        eukulele_dir = os.path.join(OUTPUTDIR, "CAG_eukulele_pleuro")
    #conda:
    #    os.path.join("..","envs","eukulele-env.yaml")
    shell:
        """
        mkdir -p {params.eukulele_dir}
        EUKulele --mets_or_mags mets --sample_dir {params.sampledir} --p_ext ".pep" --reference_dir ../2021-09-ALOHA/pleuromamma_marmmetsp_eukulele -o {params.eukulele_dir}
        touch {output.eukulele_done}
        """

rule eukulele_cag:
    input:
        pep = (lambda filename: expand(os.path.join(OUTPUTDIR, "intermediate-files",
                         "04-compare",\
                           "08-CAG-proteins", "{assembly}_CAG.fasta.transdecoder.pep"),
               assembly = assemblygroups)),
        salmon = (lambda filename: expand(os.path.join(OUTPUTDIR, "intermediate-files",\
                 "04-compare",
             "09-CAG-mapping",\
                     "salmon", "{assembly}_quant", "quant.sf"),
               assembly = assemblygroups))
    output:
        eukulele_directory = directory(os.path.join(OUTPUTDIR, "CAG_eukulele")),
        eukulele_done = os.path.join(OUTPUTDIR, "CAG_eukulele", "EUKulele_done.txt")
        #tax_est = os.path.join(OUTPUTDIR, "eukulele_CAG", "taxonomy_estimation", "{CAG}_)
    params:
        sampledir = os.path.join(OUTPUTDIR, "intermediate-files",
                                             "04-compare",\
                           "08-CAG-proteins")
    #conda:
    #    os.path.join("..","envs","eukulele-env.yaml")
    shell:
        """
        EUKulele --mets_or_mags mets --sample_dir {params.sampledir} --p_ext ".pep" --reference_dir ~/marmmetsp
    touch {output.eukulele_done}
        """
