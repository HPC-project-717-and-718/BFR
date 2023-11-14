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

typedef struct{
    int index_of_cset_1;
    int index_of_cset_2;
    double distance;
}hierc_element;

typedef struct {
    hierc_element* data;
    int size;
    int capacity;
} PriorityQueue;

RetainedSet init_retained_set();
CompressedSets init_compressed_sets();
Cluster * init_cluster(int k);
PriorityQueue createPriorityQueue(int capacity);


void print_clusters(Cluster * clusters);
void print_compressedsets(CompressedSets C);
void print_retainedset(RetainedSet R);


void add_point_to_retained_set(RetainedSet * R, Point p);
void update_cluster(Cluster * cluster, Point p);
void update_centroids(Cluster ** clusters, int number_of_clusters);
void merge_compressedsets_and_miniclusters(CompressedSets * C, Cluster * miniclusters, int number_of_miniclusters);
data_streamer data_streamer_Init(char * file_name, char * mode);
// void add_compressedset(CompressedSets * C, CompressedSet C1, bool * cset_validity);
// void remove_compressedset(CompressedSets * C, int i, int j, bool * cset_validity);
CompressedSet merge_compressedsets(CompressedSet C1, CompressedSet C2);
bool tightness_evaluation_cset(CompressedSet C);

bool add_to_pqueue(PriorityQueue * pq, hierc_element element);
bool pop_from_pqueue(PriorityQueue * pq);
hierc_element * top_from_pqueue(PriorityQueue pq);
// bool remove_from_pqueue(PriorityQueue* pq, int index1, int index2);
#endif