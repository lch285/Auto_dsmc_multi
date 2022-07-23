#PBS -S /bin/bash
#PBS -l walltime=96:00:00
#PBS -l select=1:ncpus=40:mpiprocs=40:model=sky_ele
#PBS -N Automated
#PBS -j oe
#PBS -q long
#PBS -m a
#PBS -m b
#PBS -m e


#module avail

module purge
module load gcc/8.4
module load comp-intel/2020.4.304
module load mpi-intel/2020.0.166
module load python3/3.9.5
module list

export jobid=`echo $PBS_JOBID | awk -F . '{print $1}'`

python settings.py > Automated.$jobid
