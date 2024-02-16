//
// Created by francesco.boldrin on 2/7/24.
//
#ifndef BFR_STRUCTURES
#define BFR_STRUCTURES

#input <stdio.h>
#input <stdlib.h>
#input <mpi.h>

#define M 3  // Assuming M is the size of the 'coords' array
// TODO: import the number of dimension M

typedef struct {
    double coords[M];
    int cluster;
} Point;