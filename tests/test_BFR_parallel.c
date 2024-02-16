#include <stdio.h>
#include "../src/BFR_parallel.h"
#include <time.h>

void test_hierchical_clustering_multi() {
    // Create a sample CompressedSets object
    CompressedSets C;
    // Initialize C with some values
    for (int i = 0; i < 10; i++) {
        C.sets[i].size = 10;
        for (int j = 0; j < 10; j++) {
            C.sets[i].elements[j] = i * 10 + j;
        }
    }
    
    // Call the hierchical_clustering_multi function
    hierchical_clustering_multi(&C);
    
    // Add assertions to verify the correctness of the function
    

    // print the results
    for (int i = 0; i < 10; i++) {
        printf("Cluster %d: ", i);
        for (int j = 0; j < C.sets[i].size; j++) {
            printf("%d ", C.sets[i].elements[j]);
        }
        printf("\n");
    }
}

void test_initClustersWithCentroids(){
    // create a set of clusters K
    int K = 10;
    int n = 1000;
    int d = 10;

    Cluster *clusters = (Cluster *)malloc(K * sizeof(Cluster));

    // create a list of random points in d dimensions
    Point *points = (Point *)malloc(n * sizeof(Point));

    for (int i = 0; i < n; i++) {
        points[i].coordinates = (double *)malloc(d * sizeof(double));
        for (int j = 0; j < d; j++) {
            srand(time(NULL));
            points[i].coordinates[j] = ((double) rand() %%1000)/ 100;
        }
    }


    // call the initClustersWithCentroids function
    initClustersWithCentroids(points, n, d, K, clusters);
    

    // Add assertions to verify the correctness of the function

    // print the results
    for (int i = 0; i < K; i++) {
        printf("Cluster %d: ", i);
        for (int j = 0; j < clusters[i].size; j++) {
            printf("%d ", clusters[i].elements[j]);
        }
        printf("\n");
    }
}


void 

int main() {


    return 0;
}