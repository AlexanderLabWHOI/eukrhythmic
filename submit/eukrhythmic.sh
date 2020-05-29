#!/bin/bash
#SBATCH partition=compute
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=24:00:00

alias eukrhythmic='./bin/eukrhythmic.sh'
eukrhythmic --use-sample --slurm
