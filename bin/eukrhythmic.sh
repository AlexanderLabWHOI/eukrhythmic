#!/bin/bash
#SBATCH --qos=unlim
#SBATCH --time=5000
#SBATCH --partition=compute
#SBATCH --mem=100gb

## PARAMETERS THAT MAY BE CHANGED FROM THE COMMAND LINE ##
SUBROUTINE="all" # --subroutine
JOBS=500 # --jobs; number of simultaneous Snakemake jobs
JOBNAMEARG="" # -n / --job-name
METATSAMPLEARG="" # -s / --sample-file-name
OUTDIRARG="" # -o / --out-dir
INDIRARG="" # -i / --in-dir
SCRATCHDIRARG="" # -c / --scratch-dir
CHECKQUALFLAG=1 # -q / --check-quality (Boolean)
RUNBBMAPFLAG=0 # -b / --run-bbmap (Boolean; set at the same time as spikefile)
RUNBBMAPARG="" # -b / --run-bbmap (Boolean; set at the same time as spikefile)
SLURMFLAG=0 # -l / --slurm (Boolean)
GENFILEFLAG=0 # -g / --generate-file (Boolean)
USESAMPFLAG=0 # --use-sample (Boolean)
DRYRUNFLAG=0 # --dry-run (Boolean)

# If they are not present in the config file presently, they will be empty strings.
if [[ "$(cat config.yaml | grep jobname | cut -d ":" -f 2 | cut -d " " -f 2 )" != "" ]]; then export JOBNAME="$(cat config.yaml | grep jobname | cut -d ":" -f 2 | cut -d " " -f 2)"; fi 2> err.log
if [[ "$(cat config.yaml | grep outputDIR | cut -d ":" -f 2 | cut -d " " -f 2)" != "" ]]; then export OUTDIRARG="$(cat config.yaml | grep outputDIR | cut -d ":" -f 2 | cut -d " " -f 2)"; fi 2> err.log
if [[ "$(cat config.yaml | grep scratch | cut -d ":" -f 2 | cut -d " " -f 2)" != "" ]]; then export SCRATCHDIRARG="$(cat config.yaml | grep scratch | cut -d ":" -f 2 | cut -d " " -f 2)"; fi 2> err.log
if [[ "$(cat config.yaml | grep metaT_sample | cut -d ":" -f 2 | cut -d " " -f 2)" != "" ]]; then export METATSAMPLEARG="$(cat config.yaml | grep metaT_sample | cut -d ":" -f 2 | cut -d " " -f 2)"; fi 2> err.log
if [[ "$(cat config.yaml | grep inputDIR | cut -d ":" -f 2 | cut -d " " -f 2)" != "" ]]; then export INDIRARG="$(cat config.yaml | grep inputDIR | cut -d ":" -f 2 | cut -d " " -f 2)"; fi 2> err.log
if [[ "$(cat config.yaml | grep checkqual | cut -d ":" -f 2 | cut -d " " -f 2)" != "" ]]; then export CHECKQUALFLAG="$(cat config.yaml | grep checkqual | cut -d ":" -f 2 | cut -d " " -f 2)"; fi 2> err.log
if [[ "$(cat config.yaml | grep runbbmap | cut -d ":" -f 2 | cut -d " " -f 2)" != "" ]]; then export CHECKQUALFLAG="$(cat config.yaml | grep runbbmap | cut -d ":" -f 2 | cut -d " " -f 2)"; fi 2> err.log
if [[ "$(cat config.yaml | grep spikefile | cut -d ":" -f 2 | cut -d " " -f 2)" != "" ]]; then export RUNBBMAPARG="$(cat config.yaml | grep spikefile | cut -d ":" -f 2 | cut -d " " -f 2)"; fi 2> err.log

# Baseline code for command parser borrowed from https://medium.com/@Drew_Stokes/bash-argument-parsing-54f3b81a6a8f
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -j|--jobs)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        JOBS=$2
        shift 2
      else
        echo "Error: No number of jobs ($1) provided." >&2
        exit 1
      fi
      ;;
    -t|--subroutine)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        SUBROUTINE=$2
        shift 2
      else
        echo "Error: No subroutine ($2) provided." >&2
        exit 1
      fi
      ;;
    -q|--check-quality)
      CHECKQUALFLAG=1
      sed -i "/checkqual/c\checkqual: $CHECKQUALFLAG" config.yaml
      shift
      ;;
    -g|--generate-file)
      GENFILEFLAG=1
      shift
      ;;
    --use-sample)
      USESAMPFLAG=1
      shift
      ;;
    -np|--dry-run)
      DRYRUNFLAG=1
      shift
      ;;
    -l|--slurm)
      SLURMFLAG=1
      shift
      ;;
    -n|--job-name)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        JOBNAMEARG=$2
        sed -i "/jobname/c\jobname: $JOBNAMEARG" config.yaml
        shift 2
      else
        echo "Error: No job name ($1) provided." >&2
        exit 1
      fi
      ;;
    -s|--sample-file-name)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        METATSAMPLEARG=$2
        sed -i "/metaT_sample/c\metaT_sample: $METATSAMPLEARG" config.yaml
        shift 2
      else
        echo "Error: No sample file name ($1) provided." >&2
        exit 1
      fi
      ;;
    -o|--out-dir)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        OUTDIRARG=$2
        sed -i "/outputDIR/c\outputDIR: $OUTDIRARG" config.yaml
        shift 2
      else
        echo "Error: No output directory ($1) specified." >&2
        exit 1
      fi
      ;;
    -i|--in-dir)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        INDIRARG=$2
        sed -i "/inputDIR/c\inputDIR: $INDIRARG" config.yaml
        shift 2
      else
        echo "Error: No input directory ($1) specified." >&2
        exit 1
      fi
      ;;
    -c|--scratch-dir)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        SCRATCHDIRARG=$2
        sed -i "/scratch/c\scratch: $SCRATCHDIRARG" config.yaml
        shift 2
      else
        echo "Error: No input directory ($1) specified." >&2
        exit 1
      fi
      ;;
    -b|--run-bbmap)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        RUNBBMAPARG=$2
        sed -i "/spikefile/c\spikefile: $RUNBBMAPARG" config.yaml
        RUNBBMAPFLAG=1
        sed -i "/runbbmap/c\runbbmap: $RUNBBMAPFLAG" config.yaml
        shift 2
      else
        echo "Error: You specified to run bbmap ($1), but did not list a spike file." >&2
        exit 1
      fi
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done

if [ METATSAMPLEARG == "" ] || [[ GENFILEFLAG -eq 1 ]]
then
    echo "Generating sample file..."
    python scripts/autogenerate_metaT_sample.py metasample.txt > /dev/null 2>&1
    sed -i "/metaT_sample/c\metaT_sample: input/metasample.txt" config.yaml
fi

if [[ USESAMPFLAG -eq 1 ]]
then
    echo "Using sample data and configuration..."
    cp static/config.yaml .
fi

python scripts/importworkspace.py

if [[ DRYRUNFLAG -eq 1 ]] && [[ SLURMFLAG -ne 1 ]]
then
    echo "Running a dry run."
    snakemake -np -s eukrhythmic --rerun-incomplete --use-conda -c 1
else
    if [[ SLURMFLAG -eq 1 ]]
    then
        echo "Running on SLURM."
        if [[ DRYRUNFLAG -eq 1 ]]
        then
            submit/eukrhythmic $SUBROUTINE --jobs $JOBS -np
        else
            submit/eukrhythmic $SUBROUTINE --jobs $JOBS
        fi
    else
        echo "Running locally."
        snakemake -s eukrhythmic --rerun-incomplete --jobs 100 --use-conda
    fi
fi
