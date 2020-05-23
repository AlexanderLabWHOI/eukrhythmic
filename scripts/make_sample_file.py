import os
import sys
import yaml
import pandas as pd 
import numpy as np

with open("config.yaml", "r") as configfile:
    config = yaml.load(configfile)

INPUTDIR = config["inputDIR"]

## This is assuming you already have an index file
indexfile = config["metaT_sample"] # the current index file
indexfile = pd.read_csv(indexfile, sep = "\t")
indexfile["FastqFile"] = [""] * len(indexfile.index)
sampleIDs = indexfile.SampleID
samplefiles = os.listdir(INPUTDIR)

extension = ".fastq.gz" # the extension our sample fastq files have
inputargs = sys.argv
reverse_1 = "_1" # indicates we're telling the difference between forward and reverse
reverse_2 = "_2" # indicates we're telling the difference between forward and reverse
reverse_3 = "_R" # indicates we're telling the difference between forward and reverse

if len(inputargs) > 0:
    reverse_1 = inputargs[0]
if len(inputargs) > 1:
    reverse_2 = inputargs[1]

for tt in range(len(sampleIDs)):
    # get the fastq name if this is one of the files we're looking for
    fastq_name = [g.split(extension)[0] for g in samplefiles if sampleIDs[tt] in g]
    if len(fastq_name) != 0:
        indexfile["FastqFile"][tt] = fastq_name[0].split(reverse_1)[0]
        if reverse_2 in fastq_name[0]:
            indexfile["FastqFile"][tt] = fastq_name[0].split(reverse_2)[0]
        elif reverse_3 in fastq_name[0]:
            indexfile["FastqFile"][tt] = fastq_name[0].split(reverse_3)[0]
    else:
        print(sampleIDs[tt])    

indexfile = indexfile[indexfile.SampleID != "SH265"]    
indexfile = indexfile[indexfile.FastqFile != ""]    
indexfile.to_csv(path_or_buf = config["metaT_sample"], sep = "\t")
