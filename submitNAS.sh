#PBS -S /bin/bash
#PBS -l walltime=120:00:00
#PBS -l select=nodes:ncpus=processors:mpiprocs=processors:model=sky_ele
#PBS -N slurmtemp
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

module list

mpiexec -np totalprocessors /u/lchacon/git/dsmc/src/spa_mpi -in dsmc.input
