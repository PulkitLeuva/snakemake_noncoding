#!/bin/bash 
#SBATCH --job-name=non_coding #Job name to be displayed by for example squeue 
#SBATCH --output=non_coding%j.out #Path to the file where the job (error) output is written to 
#SBATCH --nodes=1 #Number of nodes. Multiple nodes are only useful for jobs with distributed memory (e.g. MPI). 
#SBATCH --ntasks-per-node=20 #Number of (MPI) processes per node. More than one useful only  for MPI jobs. Maximum number depends nodes (number of cores) 
#SBATCH --mem=300GB  #Memory (RAM) per node. Number followed by unit prefix. 
#SBATCH --partition=newcompute #Partition/queue in which run the job. 

cd /lustre/pulkit.h/snakemake_local/non_coding_rna_snakemake
snakemake --unlock
snakemake --cores all --forceall