import pandas as pd
import numpy as np
import os
import sys
import yaml 

with open("config.yaml", "r") as configfile:
    config = yaml.load(configfile)

indexfile = config["metaT_sample"]
samplenames = list(pd.read_csv(indexfile, sep = "\t").SampleID)

# The purpose of this script is to create a tab-delimited
# file containing the names of the data sources from
# the input table and assigning each to an 
# assembly group. 

# Later, this will be adapted to take into account
# sourmash comparisons (I don't yet know how to 
# incorporate this into the snakemake workflow, but
# I will get there!). 

assemblygroups = []
assemblygroups.append(samplenames) # this will be better later; for now every sample is in assembly group 1

assemblyfile = pd.DataFrame({'AssemblyGroup': [], \
                             'SampleName': []})

for a in range(0,len(assemblygroups)):
    groupname = str(a+1)
    thisgroup = pd.DataFrame({'AssemblyGroup': [groupname] * len(assemblygroups[a]), \
                              'SampleName': assemblygroups[a]})
    assemblyfile = assemblyfile.append(thisgroup)

assemblyfile.to_csv(path_or_buf = os.path.join("input", "assemblyfile.txt"), sep = "\t")
