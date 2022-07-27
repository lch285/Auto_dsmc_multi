#!/bin/bash
#SBATCH -t 6-0:00:00 # hours, min, sec
#SBATCH --job-name=automated
#SBATCH --nodes=1
#SBATCH -n 128                      #Number of cores needed for the job
#SBATCH --partition=normal    #Name of the partition. See instructions in the pdf file sent to you. 
#SBATCH --account=coa_sjpo228_uksr #Name of account to run under; If you are part of CTFL, don't change. If you are collaborator, email Savio Poovathingal first. 
#SBATCH --exclusive
echo "Job running on SLURM NODELIST: $SLURM_NODELIST "

module load ccs/conda/python/3.9.6
source ~/conda_init.sh
conda activate myenv

python automated_porosity.py
conda deactivate 
