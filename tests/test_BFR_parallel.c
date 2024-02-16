#include <stdio.h>
#include "../src/BFR_parallel.h"
#
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


void test_hierchical_clustering_par(bool DEBUG){
    int DIM = 10;
    // Create a sample CompressedSets object
    CompressedSets *C = (CompressedSets *)malloc(sizeof(CompressedSets));
    // Initialize C with some values
    C - > number_of_sets = 10;

    // create 50 compressed sets
    C->sets = (CompressedSet *)malloc(50 * sizeof(CompressedSet));

    for (int i = 0; i < 10; i++) {
        C->sets[i] = (CompressedSet)malloc(sizeof(CompressedSet));
        C->sets[i]-> number_of_points = 50;
        
        //fill sum and sum of squares
        C->sets[i]->sum = (double *)malloc(DIM * sizeof(double));
        C->sets[i]->sum_of_squares = (double *)malloc(DIM * sizeof(double));
        srand(time(NULL));
        for (int j = 0; j < 10; j++) {
            C->sets[i]->sum[j] = (rand() % 100000)/190;
            C->sets[i]->sum_of_squares[j] = (rand() % 100000)/190;
        }
    }
    
    // create a copy of the original compressed sets
    CompressedSets *C_copy = (CompressedSets *)malloc(sizeof(CompressedSets));
    C_copy->number_of_sets = C->number_of_sets;
    C_copy->sets = (CompressedSet *)malloc(C->number_of_sets * sizeof(CompressedSet));
    for (int i = 0; i < C->number_of_sets; i++) {
        C_copy->sets[i] = (CompressedSet)malloc(sizeof(CompressedSet));
        C_copy->sets[i]->number_of_points = C->sets[i]->number_of_points;
        C_copy->sets[i]->sum = (double *)malloc(DIM * sizeof(double));
        C_copy->sets[i]->sum_of_squares = (double *)malloc(DIM * sizeof(double));
        for (int j = 0; j < DIM; j++) {
            C_copy->sets[i]->sum[j] = C->sets[i]->sum[j];
            C_copy->sets[i]->sum_of_squares[j] = C->sets[i]->sum_of_squares[j];
        }
    }

    // calculate the time taken by the parallel version
    clock_t start, end;
    start = clock();
    // Call the hierchical_clustering_par function
    hierchical_clustering_par(&C);
    end = clock();

    double time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time taken by parallel version: %f\n", time_taken);

    Cluster *clusters = (Cluster *)malloc(10 * sizeof(Cluster));


    // calculate the time taken by the sequential version
    start = clock();
    // Call the hierchical_clustering function
    merge_compressedsets_and_miniclusters(&C_copy, clusters, 0);
    end = clock();

    time_taken = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time taken by sequential version: %f\n", time_taken);

    
    // Add assertions to verify the correctness of the function
    



    // print the results of the parallel version
    for (int i = 0; i < C->number_of_sets; i++) {
        printf("Compressed sets %d: ", i);
        for (int j = 0; j < C->sets[i]->number_of_points; j++) {
            printf("%d ", C->sets[i]->points[j]);
        }
        printf("\n");
    }


    // print the results of the sequential version
    for (int i = 0; i < C_copy->number_of_sets; i++) {
        printf("Compressed sets %d: ", i);
        for (int j = 0; j < C_copy->sets[i]->number_of_points; j++) {
            printf("%d ", C_copy->sets[i]->points[j]);
        }
        printf("\n");
    }
}

int main() {    
    test_hierchical_clustering_par(false);


    return 0;
}