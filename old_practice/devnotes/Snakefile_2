configfile: "config.yaml"

import io
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
import pathlib
from snakemake.exceptions import print_exception, WorkflowError

DATAFILE = config["metaT_sample"]
INPUTDIR1 = config["inputDIR1"]
INPUTDIR2 = config["inputDIR2"]
OUTPUTDIR = config["outputDIR"]
SCRATCHDIR = config["scratch"]
INPUTFILES = [f for f in os.listdir(INPUTDIR1) if isfile(join(INPUTDIR1, f))] + [f for f in os.listdir(INPUTDIR2) if isfile(join(INPUTDIR2, f))];

columbiasamplenames = list(pd.read_table(DATAFILE).ColumbiaID);

sample_barcodes = dict();
sample_codes = dict();

# create a dictionary that contains a list with the relevant
# information about each sample: the barcode and the L code
def get_base_names():
	for i in INPUTFILES:
		split_i = i.split("_")
		if split_i[0] in columbiasamplenames:
			sample_barcodes[split_i[0]] = split_i[1]
			sample_codes[split_i[0]] = split_i[2]

get_base_names()
print(sample_barcodes.values())

include: "modules/fastqc-snake"
include: "modules/trimmomatic-snake"
include: "modules/trinity-snake"

rule all:
	input:
		# TRIMMOMATIC OUTPUTS
		trimmed = expand("{base}/firstrim/{columbia}_{barcode}_{code}_{num}.trimmed.fastq.gz", base = OUTPUTDIR, columbia = sample_barcodes.keys(), barcode = sample_barcodes.values(), code = sample_codes.values(), num = [1,2])
 
