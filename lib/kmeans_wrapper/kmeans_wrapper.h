#ifndef KMEANS_WRAPPER
#define KMEANS_WRAPPER


#include "../bfr_structures/bfr_structures.h"
#include "../kmeans/kmeans.h"


double mahalanobis_distance(Cluster c, Point p);
double distance(const Pointer a, const Pointer b);


Cluster * cluster_retained_set(RetainedSet * R, int *k);


#endif