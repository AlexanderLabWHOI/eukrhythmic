import os
import sys
import yaml
import pandas as pd 
import numpy as np

with open("config.yaml", "r") as configfile:
    config = yaml.load(configfile)

INPUTDIR = config["inputDIR"]

inputargs = sys.argv

## CHECK THAT THE INPUT IS CORRECT ##
if (len(inputargs) < 2):
    print("You need to provide the name of the sample file you wish to generate \
           as an argument. You may also provide the extension as a second argument;\
           the default is .fastq.gz. You can provide a special indication of suffix\
           or forward/reverse reads as a third and fourth argument. The default is _1/_2, or\
           the full name if no specification is made. A fifth argument may indicate an additional\
           suffix, such as 001.")
    sys.exit()
else:
    print("Ready to go.")
    print(inputargs)
    
outfilename = os.path.join("input", inputargs[1].split("/")[0])
if outfilename.split(".")[len(outfilename.split(".")) - 1] != "txt":
    outfilename = outfilename + ".txt"
candidatefiles = os.listdir(INPUTDIR)

## SPECIFY FILE EXTENSION ##
extension = ".fastq.gz"
if len(inputargs) > 2:
    extension = inputargs[2]
# filter directory by files containing proper extension.
sequencefiles = [curr for curr in candidatefiles if extension in curr]
print([curr.split("_")[0] for curr in sequencefiles])

## SPECIFY FORWARD AND REVERSE SPECIFICATION ##
forwardspec = "_1"
reversespec = "_2"
if len(inputargs) > 4:
    forwardspec = inputargs[3]
    reversespec = inputargs[4]
    
## SPECIFY ADDITIONAL SPLIT ##
singlespec = extension
if len(inputargs) > 5:
    singlespec = inputargs[5]
    
## SPECIFY FILE NAMES ##
fastqfiles_forward = set([curr.split(forwardspec)[0] for curr in sequencefiles if forwardspec in curr])
fastqfiles_reverse = set([curr.split(reversespec)[0] for curr in sequencefiles if reversespec in curr])
fastqfiles_single = set([curr.split(singlespec)[0] for curr in sequencefiles if (forwardspec not in curr) & (reversespec not in curr)])

## COMBINE ALL FILENAMES ##
fastalist = list(fastqfiles_forward.intersection(fastqfiles_reverse)) + list(fastqfiles_single)
print(fastalist)

## ATTEMPT TO MAKE SHORTER SAMPLE ID NAME ##
samplelist = [curr.split("_")[0] for curr in fastalist] # try taking first token from name
if len(set(samplelist)) != len(set(fastalist)):
    for s in range(len(samplelist)):
        if samplelist.count(samplelist[s]) > 1:
            samplelist[s] = str(samplelist[s]) + "_" + str(samplelist.count(samplelist[s]))
    if len(set(samplelist)) != len(set(fastalist)):        
        samplelist = fastalist # otherwise just set sample list to fasta names
    
metaT = pd.DataFrame({"SampleName": fastalist, \
              "SampleID": samplelist, \
              "AssemblyGroup": samplelist, \
              "FastqFile": fastalist,})
metaT = metaT.sort_values(by="SampleName")

metaT.to_csv(path_or_buf = outfilename, index=False, sep = "\t")
