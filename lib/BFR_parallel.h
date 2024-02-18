#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../lib/bfr_helper_functions/bfr_helper_functions.h"


// TODO: Add the necessary declarations for the parallel version of the BFR algorithm

void hierachical_clustering_thr(CompressedSets * C);
bool hierachical_clust_parallel(CompressedSets * C, int rank, int size);
bool UpdateCentroidsMultithr(Cluster *clusters, int index, int number_of_clusters, int dimension);
void add_cluster_to_compressed_sets(CompressedSets *compressedSets, Cluster c);
void add_point_to_compressed_sets(CompressedSet *compressedSets, Point p);
void add_point_to_cluster(Cluster *c, Point p);
void remove_cset_from_compressed_sets(CompressedSets *C, CompressedSet * c, int index);

bool primary_compression_criteria(Cluster *clusters, Cluster *clusters_copy, Point p) ;
void secondary_compression_criteria(Cluster *clusters, RetainedSet *retainedSet, CompressedSets *compressedSets, int rank, int size, MPI_Datatype PointType);

CompressedSet *merge_cset(CompressedSet *c1, CompressedSet *c2);

Cluster *initClustersWithCentroids(Point *data_buffer, int size, int k, int dimension);


