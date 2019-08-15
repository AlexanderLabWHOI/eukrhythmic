#!/bin/bash
#SBATCH --job-name=trim_test
#SBATCH --output=whathappened.txt
#SBATCH --ntasks=1
#SBATCH --time=20:00
#SBATCH --mem-per-cpu=100000

snakemake   \
        /vortexfs1/omics/alexander/data/WAP/arianna-snakemake-output/firsttrim/SH410_CTTGTA_L005_1.trimmed.fastq.gz --use-conda -s ../modules/trimmomatic-snake 
