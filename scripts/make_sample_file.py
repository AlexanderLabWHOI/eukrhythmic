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
#print(sampleIDs)

test = False # get rid of this later; special case in test data
if test:
    sampleIDs = [s + "_test" for s in sampleIDs]

extension = ".fastq.gz" # the extension our sample fastq files have
reverse = "_R" # indicates we're telling the difference between forward and reverse

for tt in range(len(sampleIDs)):
    # get the fastq name if this is one of the files we're looking for
    fastq_name = [g.split(extension)[0] for g in samplefiles if sampleIDs[tt] in g]
    if len(fastq_name) != 0:
        indexfile["FastqFile"][tt] = fastq_name[0].split(reverse)[0]
    else:
        print(sampleIDs[tt])    

indexfile = indexfile[indexfile.SampleID != "SH265"]    
indexfile = indexfile[indexfile.FastqFile != ""]    
indexfile.to_csv(path_or_buf = config["metaT_sample"], sep = "\t")
