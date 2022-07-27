#!/bin/bash
#SBATCH -t 5-0:00:00 # hours, min, sec
#SBATCH --job-name=np_101
#SBATCH --nodes=nodes
#SBATCH -n totalprocessors                      #Number of cores needed for the job
#SBATCH --partition=normal    #Name of the partition. See instructions in the pdf file sent to you.
#SBATCH --account=coa_sjpo228_uksr #Name of account to run under; If you are part of CTFL, don't change. If you are collaborator, email Savio Poovathingal first.
#SBATCH --exclusive
echo "Job running on SLURM NODELIST: $SLURM_NODELIST "

mpirun -np totalprocessors /project/sjpo228_uksr/BenDeaton/git/dsmc/src/spa_mpi -in dsmc.input
