#!/bin/bash
#SBATCH --qos=unlim

snakemake   \
    --jobs 100 --rerun-incomplete --use-conda --rerun-incomplete --cluster-config cluster.yaml --cluster "sbatch --parsable --qos=unlim --partition={cluster.queue} --job-name=singlecell.{rule}.{wildcards} --mem={cluster.mem}gb --time={cluster.time} --ntasks={cluster.threads} --nodes={cluster.nodes}"


