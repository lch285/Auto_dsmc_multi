#!/bin/bash
#SBATCH -t 03-00:00:00 # hours, min, sec
#SBATCH --job-name=np_101
#SBATCH -N nodes
#SBATCH --partition=CAL48M192_L    #Name of the partition. See instructions in the pdf file sent to you.
#SBATCH --account=col_sjpo228_uksr #Name of account to run under; If you are part of CTFL, don't change. If you are collaborator, email Savio Poovathingal first.
#SBATCH --exclusive
echo "Job running on SLURM NODELIST: $SLURM_NODELIST "
​
mpirun -np totalprocessors /project/sjpo228_uksr/LuisChacon/git/dsmc/src/spa_mpi -in dsmc.input