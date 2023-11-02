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

Cluster * init_cluster(){
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
    Cluster * clusters = malloc(K * sizeof(Cluster));
    
    if (clusters == NULL){
        printf("Error: could not allocate memory\n");
        exit(1);
    }

    int i = 0;
    for (i = 0; i < K; i++){
        Cluster  c;
        int j = 0;
        for (j = 0; j < M; j++){
            c.centroid.coords[j] = 0.;
            c.sum[j] = 0;
            c.sum_squares[j] = 0;
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
            printf("%d ", clusters[i].sum[j]);
        }

        printf("\n");

        printf("Cluster %d sum_squares: ", i);

        for (j = 0; j < M; j++){
            printf("%d ", clusters[i].sum_squares[j]);
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
