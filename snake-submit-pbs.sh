snakemake --rerun-incomplete --jobs 50 --use-conda --cluster-config cluster-pbs.yaml --cluster "qsub -N {rule}.{wildcards} -l vmem={cluster.mem}gb,walltime={cluster.time},nodes={cluster.nodes}:ppn={cluster.cpupertask}"

