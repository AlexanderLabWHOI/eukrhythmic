#!/bin/bash -l        
#PBS -l walltime=48:00:00,nodes=1:ppn=1,mem=5gb 

jobname=$(cat config.yaml | grep jobname | cut -d ":" -f 2 | cut -d " " -f 2)
outdir=$(cat config.yaml | grep outputDIR | cut -d ":" -f 2 | cut -d " " -f 2)

SNAKEFILE=$1 # to run all of eukrhythmic, this is "eukrhythmic"
JOBS=$2 # should default to 500
RERUN_INCOMPLETE=$3

snakemake -s $SNAKEFILE $RERUN_INCOMPLETE --jobs $JOBS --use-conda --cluster-config cluster-pbs.yaml --cluster 
"qsub -N $jobname.{rule}.{wildcards} -l vmem={cluster.mem}gb,walltime={cluster.time},nodes={cluster.nodes}:ppn={cluster.cpupertask}"
