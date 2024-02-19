#!/bin/bash

#PBS -l select=5:ncpus=4:mem=2gb

# set max execution time
#PBS -l walltime=0:02:00 -l place=pack:excl

# imposta la coda di esecuzione
#PBS -q short_cpuQ
module load mpich-3.2
module load valgrind-3.15.0 
# OMP_NUM_THREADS=4 mpirun.actual -n 5 valgrind --track-origins=yes --leak-check=full --verbose HPC/BFR/src/BFR_parallel HPC/BFR/data/synthetic/synthetic_d3_6000points.txt > HPC/BFR/output.txt 2>&1 # Full debug setup
# OMP_NUM_THREADS=4 mpirun.actual -n 5 valgrind HPC/BFR/src/BFR_parallel HPC/BFR/data/synthetic/synthetic_d3_6000points.txt > HPC/BFR/output.txt 2>&1 # Debug setup
OMP_NUM_THREADS=4 mpirun.actual -n 5 HPC/BFR/src/BFR_parallel HPC/BFR/data/synthetic/synthetic_d2_4000points.txt 

# poi esegui "qsub bfr_parallel.sh"
# non eseguire MAI manualmente!!
