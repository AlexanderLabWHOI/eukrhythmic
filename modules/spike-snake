configfile: "config.yaml"

import io
import os
import pathlib
import pandas as pd
from snakemake.exceptions import print_exception, WorkflowError
import sys
sys.path.insert(1, '../scripts')
from importworkspace import *

rule salmon_rename:
    input:
        input = (lambda filename: expand(os.path.join(OUTPUTDIR, "salmon_{folder}", "against_mega", 
                                    "salmon_quant_assembly_{assembly}", "quant.sf"), 
                       assembly = filename.assembly, folder = filename.folder))
    output:
        out = os.path.join(OUTPUTDIR, "salmon_{folder}", "against_mega", 
                           "salmon_quant_assembly_{assembly}", "quant.sf_cleaned")
    params:
        extra = ""
    shell:
        """
        for i in {input.input}; do awk -F, -v OFS=, 'NR==1{{split(FILENAME,a,"/quant.sf");$2= a[1] ""}}1' ${{i}} | awk '{{gsub(".*salmon_quant_assembly_","",$5)}}1' > {output.out}; done
        """

rule salmon_table:
    input:
        input = (lambda filename: expand(os.path.join(OUTPUTDIR, "salmon_{folder}", "against_mega", 
                                                      "salmon_quant_assembly_{assembly}", "quant.sf_cleaned"), 
                                         assembly = assemblygroups, folder = filename.folder))
    output:
        out = os.path.join(OUTPUTDIR, "salmon_{folder}", "mergedtable.tab")
    shell:
        """
        awk '
            {{samples[$1] = samples[$1] OFS $NF}}
            END {{
                print "Name", samples["Name"]
                delete samples["Name"]
                for (name in samples) print name, samples[name]
            }}
        ' {input.input} > {output.out}
        """

rule copies:
    input:
        countfile = os.path.join(OUTPUTDIR, "salmon_{folder}", "mergedtable.tab"),
        copiestable = SPIKETABLE
    output:
        out = os.path.join(OUTPUTDIR, "salmon_{folder}", "copiesperL.tab")
    params:
        spikecopies = 500000,
        IDs = assemblygroups
    run:
        final = pd.DataFrame()
        spike=pd.read_csv(input.copiestable, sep=' ')
        counts=pd.read_csv(input.countfile, sep=' ')
        for i in params.IDs:
            df = (spike.loc[spike['ID'] == (i)])
            #Spike calculation reference: https://www.protocols.io/view/internal-genomic-dna-standard-for-quantitative-met-jxdcpi6?step=12
            calc = (counts[i] * ((float((df['TotalReads'])) * 
                                  (params.spikecopies / (float((df['bbmap-total']))))) / 
                                 (float((df['TotalReads']))))) / float((df['VolumeFiltered']))
            colNames = calc.name
            rowNames = counts.iloc[:,0]
            calc.index = rowNames
            calc.column_name = colNames
            final = final.append(calc)
        final_transpose = final.T
        print(final_transpose)
        final_transpose.to_csv(output.out, sep = ' ', index=True)
