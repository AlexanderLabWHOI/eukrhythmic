#!/bin/bash
#SBATCH --qos=unlim
#SBATCH --mem=100gb
#SBATCH --partition=compute

snakemake   \
    --unlock --jobs 100 --rerun-incomplete --use-conda --rerun-incomplete --cluster-config cluster.yaml --cluster "sbatch --parsable --qos=unlim --partition={cluster.queue} --job-name=batseukrhythmic.{rule}.{wildcards} --mem={cluster.mem}gb --time={cluster.time} --ntasks={cluster.threads} --nodes={cluster.nodes}"
