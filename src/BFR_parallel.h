#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include <mpi.h>
#include <omp.h>
#include "../lib/kmeans_wrapper/kmeans_wrapper.h"

// TODO: Add the necessary declarations for the parallel version of the BFR algorithm

void hierachical_clustering_thr(CompressedSets * C);
bool UpdateCentroidsMultithr(Cluster * c);
void add_cluster_to_compressed_sets(CompressedSet *compressedSets, Cluster c);
void add_point_to_compressed_sets(CompressedSet *compressedSets, Point p);
void add_point_to_cluster(Cluster *c, Point p);
void remove_cset_from_compressed_sets(CompressedSet *C, CompressedSet * c);

bool primary_compression_criteria(Cluster *clusters, Point p) ;
bool secondary_compression_criteria(Cluster *clusters, CompressedSets *compressedSets, RetainedSet *retainedSet);

CompressedSet *merge_cset(CompressedSet *c1, CompressedSet *c2);

Cluster *initClustersWithCentroids(Point *data_buffer, int size, int K, int DIMENSION);


