#ifndef BFR_STRUCTURES
#define BFR_STRUCTURES


#include "../definitions.h"


typedef struct {
    double coords[M];
    int cluster;
} Point;

typedef struct {
    Point centroid;
    int size;
    double sum[M];
    double sum_squares[M];
    int index;
} Cluster;

// TODO: consider making the points field a linked list
typedef struct {
    Point * points;
    int number_of_points;
} RetainedSet;

typedef struct {
    int number_of_points;
    int sum[M];
    int sum_square[M];
} CompressedSet;

typedef struct {
    CompressedSet * sets;
    int number_of_sets;
} CompressedSets;


RetainedSet init_retained_set();
CompressedSets init_compressed_sets();
Cluster * init_cluster(int k);


void print_clusters(Cluster * clusters);
void print_compressedsets(CompressedSets C);
void print_retainedset(RetainedSet R);


void add_point_to_retained_set(RetainedSet * R, Point p);
void update_cluster(Cluster * cluster, Point p);
void update_centroids(Cluster ** clusters, int number_of_clusters);
void merge_compressedsets_and_miniclusters(CompressedSets * C, Cluster * miniclusters, int number_of_miniclusters);

#endif