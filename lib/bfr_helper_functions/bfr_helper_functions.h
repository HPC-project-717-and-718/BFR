#ifndef KMEANS_WRAPPER
#define KMEANS_WRAPPER


#include "../bfr_structures/bfr_structures.h"


double mahalanobis_distance(Cluster c, Point p);
double distance(const Pointer a, const Pointer b);


Cluster * cluster_retained_set(RetainedSet * R, int *k);
kmeans_config init_kmeans_config(int k, RetainedSet * R, bool parallel, int rank, int size);

bool load_data_buffer(data_streamer stream_cursor, Point * data_buffer, int max_size_of_buffer, long * size_of_data_buffer);

#endif