#!/bin/bash

#PBS -l select=1:ncpus=1:mem=2gb

# set max execution time
#PBS -l walltime=0:02:00 -l place=pack:excl

# imposta la coda di esecuzione
#PBS -q short_cpuQ
module load mpich-3.2
mpirun.actual -n 1 HPC/BFR/src/bfr_serial HPC/BFR/data/synthetic/synthetic_d2_4000points_2gaussians.txt

# poi esegui "qsub pingpong.sh"
# non eseguire MAI manualmente!!
