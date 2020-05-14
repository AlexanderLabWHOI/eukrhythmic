#!/bin/bash
#SBATCH --qos=unlim
#SBATCH --time=5000
#SBATCH --partition=compute
#SBATCH --mem=100gb

jobname=$(cat config.yaml | grep jobname | cut -d ":" -f 2 | cut -d " " -f 2)

snakemake  \
    --jobs 100 --rerun-incomplete --use-conda --rerun-incomplete --cluster-config cluster.yaml --cluster "sbatch --parsable --qos=unlim --partition={cluster.queue} --job-name=${jobname}.{rule}.{wildcards} --mem={cluster.mem}gb --time={cluster.time} --ntasks={cluster.threads} --nodes={cluster.nodes}"
