#!/bin/bash
#SBATCH -t 04-00:00:00 # hours, min, sec
#SBATCH --job-name=np_multi
#SBATCH -N 1                      #Number of cores needed for the job
#SBATCH --partition=CAC48M192_L    #Name of the partition. See instructions in the pdf file sent to you.
#SBATCH --account=col_sjpo228_uksr #Name of account to run under; If you are part of CTFL, don't change. If you are collaborator, email Savio Poovathingal first.
#SBATCH --exclusive
echo "Job running on SLURM NODELIST: $SLURM_NODELIST "

module load ccs/conda/python-3.7.3
conda activate myenv

#python automated_porosity.py
#python lhs_automated.py
python settings.py
conda deactivate
