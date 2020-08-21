#!/bin/bash
#SBATCH --qos=unlim
#SBATCH --time=5000
#SBATCH --partition=scavenger
#SBATCH --mem=100gb

jobname=$(cat config.yaml | grep jobname | cut -d ":" -f 2 | cut -d " " -f 2)
outdir=$(cat config.yaml | grep outputDIR | cut -d ":" -f 2 | cut -d " " -f 2)

SNAKEFILE=$1 # to run all of eukrhythmic, this is "eukrhythmic"
JOBS=$2 # should default to 500
RERUN_INCOMPLETE=$3

snakemake  \
    -s $SNAKEFILE $RERUN_INCOMPLETE --jobs $JOBS --use-conda --cluster-config cluster.yaml --cluster "sbatch --parsable --qos=unlim --partition={cluster.queue} --job-name=${jobname}.{rule}.{wildcards} --mem={cluster.mem}gb --time={cluster.time} --ntasks={cluster.tasks} --nodes={cluster.nodes} --cpus-per-task={cluster.cpupertask} --requeue"
