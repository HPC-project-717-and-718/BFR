#include "bfr_structures.h"

RetainedSet init_retained_set(){
    /*
    * Initialize retained set
    *
    * Algorithm:
    *   1. allocate memory for retained set
    *   2. initialize retained set as NULL
    *
    * Parameters:
    *   - void
    *
    * Returns:
    *   - retained set
    */
    RetainedSet R;
    R.points = NULL;
    R.number_of_points = 0;
    return R;
}

CompressedSets init_compressed_sets(){
    /*
    * Initialize compressed sets
    *
    * Algorithm:
    *   1. allocate memory for compressed sets
    *   2. initialize each compressed set as NULL
    *
    * Parameters:
    *   - void
    *
    * Returns:
    *   - array of compressed sets
    */
    CompressedSets C;
    C.sets = NULL;
    C.number_of_sets = 0;
    return C;
}

Cluster * init_cluster(int k){
    /*
    * Initialize clusters
    *
    * Algorithm:
    *   1. allocate memory for clusters
    *   2. initialize each cluster
    *
    * Parameters:
    *   - void
    *
    * Returns:
    *   - array of clusters
    */
    Cluster * clusters = malloc(k * sizeof(Cluster));
    
    if (clusters == NULL){
        printf("Error: could not allocate memory\n");
        exit(1);
    }

    int i = 0;
    for (i = 0; i < k; i++){
        Cluster  c;
        int j = 0;
        for (j = 0; j < M; j++){
            c.centroid.coords[j] = 0.;
            c.sum[j] = 0.;
            c.sum_squares[j] = 0.;
        }
        c.size = 0;
        c.index = i;
        clusters[i] = c;
    }

    return clusters;
}


void print_clusters(Cluster * clusters){
    printf("Clusters:\n");
    int i = 0;
    for (i = 0; i < K; i++){
        printf("Cluster %d: ", i);
        int j = 0;
        for (j = 0; j < M; j++){
            printf("%lf ", clusters[i].centroid.coords[j]);
        }
        printf("\n");

        printf("Cluster %d size: %d\n", i, clusters[i].size);

        printf("Cluster %d sum: ", i);

        for (j = 0; j < M; j++){
            printf("%lf ", clusters[i].sum[j]);
        }

        printf("\n");

        printf("Cluster %d sum_squares: ", i);

        for (j = 0; j < M; j++){
            printf("%lf ", clusters[i].sum_squares[j]);
        }

        printf("\n");
    }
}

void print_compressedsets(CompressedSets C){
    printf("Compressed Sets:\n");
    int i = 0;
    for(; i < C.number_of_sets; i++){
        printf("Set %d size: %d\n", i, C.sets[i].number_of_points);

        printf("Set %d sum: ", i);
        int j = 0;
        for (j = 0; j < M; j++){
            printf("%d ", C.sets[i].sum[j]);
        }
        printf("\n");

        printf("Set %d sum_squares: ", i);
        j = 0;
        for (j = 0; j < M; j++){
            printf("%d ", C.sets[i].sum_square[j]);
        }
        printf("\n");
    }
}

void print_retainedset(RetainedSet R){
    printf("Retained Set:\n");
    printf("Retained Set size: %d\n", R.number_of_points);
    int i = 0;
    for (; i < R.number_of_points; i++){
        printf("Point %d: ", i);
        int j = 0;
        for (j = 0; j < M; j++){
            printf("%lf ", R.points[i].coords[j]);
        }
        printf("\n");
    }
}


void add_point_to_retained_set(RetainedSet * R, Point p){
    (*R).number_of_points += 1;
    (*R).points = realloc(R->points, R->number_of_points * sizeof(Point));
    if (R->points == NULL){
        printf("Error: could not allocate memory\n");
        exit(1);
    }
    (*R).points[R->number_of_points - 1] = p;
}

void update_cluster(Cluster * cluster, Point p){
    /*
    * Update cluster
    *
    * Algorithm:
    *   1. update cluster size
    *   2. update cluster sum
    *   3. update cluster sum_squares
    *
    * Parameters:
    *   - cluster: cluster to update
    *   - p: point to add to cluster
    *
    * Returns:
    *   - void
    */
    cluster->size += 1;
    int i = 0;
    for (i = 0; i < M; i++){
        cluster->sum[i] += p.coords[i];
        cluster->sum_squares[i] += pow(p.coords[i], 2);
    }
}

void update_centroids(Cluster ** clusters, int number_of_clusters){
    /*
    * Update centroids of clusters
    *
    * Algorithm:
    *   1. iter over clusters
    *   2. if cluster size is greater than 1, update centroid
    *
    * Parameters:
    *   - clusters: array of clusters
    *   - number_of_clusters: number of clusters
    *
    * Returns:
    *   - void
    */

    int i;
    for (i = 0; i < number_of_clusters; i++){
        int size = (*clusters)[i].size;
        if (size >= 1){
            Point new_centroid;
            int j = 0;
            for (j = 0; j < M; j++){
                new_centroid.coords[j] = (double) (*clusters)[i].sum[j] / (double) (*clusters)[i].size; 
            }
            new_centroid.cluster = (*clusters)[i].centroid.cluster;
            (*clusters)[i].centroid = new_centroid;
        }
    }
}

void add_miniclusters_to_compressedsets(CompressedSets * C, Cluster * miniclusters, int number_of_miniclusters){
    int i;
    for (i=0; i<number_of_miniclusters; i++){
        if (miniclusters[i].size > 1){
            // add minicluster to compressed sets
            (*C).number_of_sets += 1;
            (*C).sets = realloc(C->sets, C->number_of_sets * sizeof(CompressedSet));
            if (C->sets == NULL){
                printf("Error: could not allocate memory\n");
                exit(1);
            }
            (*C).sets[C->number_of_sets - 1].number_of_points = miniclusters[i].size;
            int j;
            for (j=0; j<M; j++){
                (*C).sets[C->number_of_sets - 1].sum[j] = miniclusters[i].sum[j];
                (*C).sets[C->number_of_sets - 1].sum_square[j] = miniclusters[i].sum_squares[j];
            }
        }
    }
}

void merge_compressedsets_and_miniclusters(CompressedSets * C, Cluster * miniclusters, int number_of_miniclusters){
    add_miniclusters_to_compressedsets(C, miniclusters, number_of_miniclusters);

    // TODO: merge compressed sets using hierarchical clustering
    // priority queue pseudocode (naive implementation is also possible, I don't care tbh)
    // priority queue of compressed sets
    // while (priority queue is not empty){
    //     pop two compressed sets
    //     merge them if they fit criteria (to be discussed)
    //     push merged compressed set to priority queue, if it fits criteria
    // }
}