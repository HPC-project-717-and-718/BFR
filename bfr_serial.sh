#!/bin/bash

#PBS -l select=1:ncpus=1:mem=2gb

# set max execution time
#PBS -l walltime=0:02:00 -l place=pack:excl

# imposta la coda di esecuzione
#PBS -q short_cpuQ
module load mpich-3.2
module load valgrind-3.15.0 
# mpirun.actual -n 1 valgrind --track-origins=yes --leak-check=full --verbose HPC/BFR/src/bfr_serial HPC/BFR/data/synthetic/synthetic_d2_4000points_2gaussians.txt > HPC/BFR/output.txt 2>&1 # Full debug setup
#mpirun.actual -n 1 valgrind HPC/BFR/src/bfr_serial HPC/BFR/data/synthetic/synthetic_d2_4000points_2gaussians.txt > HPC/BFR/output.txt 2>&1 # Debug setup
mpirun.actual -n 1 time HPC/BFR/src/bfr_serial HPC/BFR/data/synthetic/synthetic_d2_4000points.txt 

# poi esegui "qsub pingpong.sh"
# non eseguire MAI manualmente!!
