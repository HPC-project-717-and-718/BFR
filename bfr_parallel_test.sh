#!/bin/bash

#PBS -l select=3:ncpus=2:mem=2gb

# set max execution time
#PBS -l walltime=0:02:00 -l place=pack:excl

# imposta la coda di esecuzione
#PBS -q short_cpuQ
module load mpich-3.2
module load valgrind-3.15.0 
# OMP_NUM_THREADS=2 mpirun.actual -n 3 valgrind --track-origins=yes --leak-check=full --verbose HPC/BFR/example2.o > HPC/BFR/output.txt 2>&1 # Full debug setup
# OMP_NUM_THREADS=2 mpirun.actual -n 3 valgrind HPC/BFR/example2.o > HPC/BFR/output.txt 2>&1 # Debug setup
OMP_NUM_THREADS=2 mpirun.actual -n 4 HPC/BFR/tests/test_BFR_parallel.o

# poi esegui "qsub bfr_parallel_test.sh"
# non eseguire MAI manualmente!!

# mpicc tests/test_BFR_parallel.c lib/bfr_structures/bfr_structures.c lib/kmeans_wrapper/kmeans_wrapper.c lib/kmeans/kmeans.c src/BFR_parallel.c -o test/test_BFR_parallel.o -lm -fopenmp