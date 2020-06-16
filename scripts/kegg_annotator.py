import pandas as pd
import sys
import os
import numpy as np
import argparse

#### PROCESSING FUNCTIONS ####
def parseFile(file_in, column_1, column_2):
    file_read = pd.read_csv(file_in, header = None, sep = "\t")
    for c in file_read.columns:
        file_read[c] = [curr.split(":")[len(curr.split(":"))-1] for curr in file_read[c]]
    return file_read

def createDictionary(list_in):
    dict_out = dict()
    for r in range(len(list_in.index)):
        curr = list_in.iloc[r,0]
        if (curr in dict_out):
            if (dict_out[curr] != None):
                dict_out[curr] = dict_out[curr].append(list_in.iloc[r,1])
            else:
                dict_out[curr] = [list_in.iloc[r,1]]
        else:
            dict_out[curr] = [list_in.iloc[r,1]]
    return dict_out

#### PROCESS KEGG PATH ARGUMENTS ####
parser = argparse.ArgumentParser(description='Read in the relevant information for your KEGG install.')
parser.add_argument('kegg_path', metavar='path', type=str,\
                    help='The path to the KEGG database. Can be specified without other paths for default paths to be inferred.')
args = parser.parse_args()

#### SPECIFIC PATH ARGUMENTS; MAY BE LEFT OUT ####
parser.add_argument('--ko-genes-list', dest='ko_genes_list', \
                    default=os.path.join(args.kegg_path, "genes", "ko", "ko_genes.list"), \
                    help='The path to the KO genes list (not required).')
parser.add_argument('--ko-module-list', dest='ko_module_list', \
                    default=os.path.join(args.kegg_path, "genes", "ko", "ko_module.list"), \
                    help='The path to the KO module list (not required).')
parser.add_argument('--ko-pathway-list', dest='ko_pathway_list', \
                    default=os.path.join(args.kegg_path, "genes", "ko", "ko_pathway.list"), \
                    help='The path to the KO pathway list (not required).')
parser.add_argument('--ko-enzyme-list', dest='ko_enzyme_list', \
                    default=os.path.join(args.kegg_path, "genes", "ko", "ko_enzyme.list"), \
                    help='The path to the KO enzyme list (not required).')
parser.add_argument('--ko', dest='ko', \
                    default=os.path.join(args.kegg_path, "genes", "ko", "ko"), \
                    help='The path to the KO file (not required).')
parser.add_argument('--module', dest='module', \
                    default=os.path.join(args.kegg_path, "module", "module"), \
                    help='The path to the module file (not required).')
parser.add_argument('--pathway', dest='pathway', \
                    default=os.path.join(args.kegg_path, "pathway", "pathway"), \
                    help='The path to the module file (not required).')

args = parser.parse_args()

#### READ IN EACH OF THE LIST FILES ####
ko_genes_list = parseFile(args.ko_genes_list, column_1 = "ko", column_2 = "genes")
ko_module_list = parseFile(args.ko_module_list, column_1 = "ko", column_2 = "module")
ko_pathway_list = parseFile(args.ko_pathway_list, column_1 = "ko", column_2 = "pathway")
ko_enzyme_list = parseFile(args.ko_enzyme_list, column_1 = "ko", column_2 = "enzyme")

#### CREATE DICTIONARY FROM EACH OF THE LIST FILES ####
ko_genes_dict = createDictionary(ko_genes_list)
ko_module_dict = createDictionary(ko_module_list)
ko_pathway_dict = createDictionary(ko_pathway_list)
ko_enzyme_dict = createDictionary(ko_enzyme_list)

#### PARSE KO FILE ####
KOfile = open(args.ko, "r")
curr_def = ""
curr_entry = ""
curr_name = ""
KO_dict = dict()
for line in KOfile:
    if "DEFINITION" in line:
        curr_def = line.split("DEFINITION")[1].strip()
    elif "ENTRY" in line:
        curr_entry = line.split("ENTRY")[1].split("KO")[0].strip()
    elif "NAME" in line:
        curr_name = line.split("NAME")[1].strip()
    elif "///" in line:
        if (not curr_def) & (not curr_entry) & (not curr_name):
            next
        elif (not not curr_entry) & (not not curr_name):
            KO_dict[curr_entry] = {"name": curr_name, "def": curr_def}
            curr_def = ""
            curr_entry = ""
            curr_name = ""
        else:
            print("Parsing error")
            
#### PARSE MODULE FILE ####
modulefile = open(args.module, "r")
curr_class = ""
curr_entry = ""
curr_name = ""
module_dict = dict()
for line in modulefile:
    if "CLASS" in line:
        curr_class = line.split("CLASS")[1].strip()
    elif "ENTRY" in line:
        curr_entry = line.split("ENTRY")[1].split("Pathway")[0].strip()
    elif "NAME" in line:
        curr_name = line.split("NAME")[1].strip()
    elif "///" in line:
        if (not curr_class) & (not curr_entry) & (not curr_name):
            next
        elif (not not curr_entry) & (not not curr_name):
            module_dict[curr_entry] = {"name": curr_name, "class": curr_class}
            curr_def = ""
            curr_entry = ""
            curr_name = ""
        else:
            print("Parsing error")

#### PARSE PATHWAY FILE ####
pathwayfile = open(args.pathway, "r")
curr_class = ""
curr_entry = ""
curr_name = ""
curr_desc = ""
pathway_dict = dict()
for line in pathwayfile:
    if "CLASS" in line:
        curr_class = line.split("CLASS")[1].strip()
    elif "ENTRY" in line:
        curr_entry = line.split("ENTRY")[1].split("Pathway")[0].strip()
    elif "DESCRIPTION" in line:
        curr_desc = line.split("DESCRIPTION")[1].strip()
    elif "NAME" in line:
        curr_name = line.split("NAME")[1].strip()
    elif "///" in line:
        if (not curr_class) & (not curr_entry) & (not curr_name):
            next
        elif (not not curr_entry) & (not not curr_name):
            pathway_dict[curr_entry] = {"name": curr_name, "desc": curr_desc, "class": curr_class}
            curr_desc = ""
            curr_entry = ""
            curr_name = ""
            cur_class = ""
        else:
            print("Parsing error")