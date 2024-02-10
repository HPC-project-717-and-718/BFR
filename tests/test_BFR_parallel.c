#include <stdio.h>
#include "../src/BFR_parallel.h"

int main() {
    // Create a sample CompressedSets object
    CompressedSets C;
    // Initialize C with some values
    for (int i = 0; i < 10; i++) {
        C.sets[i].size = 10;
        for (int j = 0; j < 10; j++) {
            C.sets[i].elements[j] = i * 10 + j;
        }
    }
    
    // Call the hierchieal_clustering_thr function
    hierchieal_clustering_thr(&C);
    
    // Add assertions to verify the correctness of the function
    

    // print the results
    for (int i = 0; i < 10; i++) {
        printf("Cluster %d: ", i);
        for (int j = 0; j < C.sets[i].size; j++) {
            printf("%d ", C.sets[i].elements[j]);
        }
        printf("\n");
    }

    return 0;
}