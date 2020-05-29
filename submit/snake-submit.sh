#!/bin/bash
#SBATCH --qos=unlim
#SBATCH --time=5000
#SBATCH --partition=scavenger
#SBATCH --mem=100gb

jobname=$(cat config.yaml | grep jobname | cut -d ":" -f 2 | cut -d " " -f 2)
outdir=$(cat config.yaml | grep outputDIR | cut -d ":" -f 2 | cut -d " " -f 2)

python scripts/writeconfig.py

snakemake  \
    --rerun-incomplete --jobs 100 --use-conda --cluster-config cluster.yaml --cluster "sbatch --parsable --qos=unlim --partition={cluster.queue} --job-name=${jobname}.{rule}.{wildcards} --mem={cluster.mem}gb --time={cluster.time} --ntasks={cluster.threads} --nodes={cluster.nodes}"
