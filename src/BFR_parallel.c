// Description: This file contains the parallel implementation of the BFR algorithm.
#include <mpi.h>
#include <omp.h>
#include "../lib/kmeans_wrapper/kmeans_wrapper.h"
#include "./BFR_parallel.h"

// TODO: there will some problem here

#define MASTER 0
#define DEBUG 1
# define K 10
# define DIMENSION 2
# define NUMBER_OF_THREADS 4

cluster *initClustersWithCentroids(Point *data_buffer, int size) {
    // TODO: implement the function, THIS IS JUST A NAIVE IMPLEMENTATION
    
    // take the first k points as centroids
    Cluster *clusters = (Cluster *)malloc(K * sizeof(Cluster));
    int i;
    for (i = 0; i < K; i++) {
        clusters[i].centroid = data_buffer[i];
        clusters[i].size = 0;
    }
    return clusters;
}

CompressedSet *merge_cset(CompressedSet c1, CompressedSet c2) {
    // allocate memory for the new compressed set
    CompressedSet * new_cset = (CompressedSet *)malloc(sizeof(CompressedSet));


    // merge the two compressed sets 
    new_cset -> number_of_points = c1 -> number_of_points + c2 -> number_of_points;
    
    # pragma omp parallel for shared(new_cset, c1, c2)
    int i;
    for ( i = 0; i < M; i++){
        new_cset -> sum[i] = c1 -> sum[i] + c2 -> sum[i];
        new_cset -> sum_square[i] = c1 -> sum_square[i] + c2 -> sum_square[i];
    }

    return new_cset;
}

void hierchieal_clustering_thr(CompressedSets * C){
    // TODO: implement this function
    // hierchieal clustering of the compressed sets

    bool stop_criteria = false;
    bool changes = false;

    float **distance_matrix =  malloc (sizeof(float *) * C->size);
    int i;
    for (i = 0; i < C->size; i++) {
        distance_matrix[i] = malloc (sizeof(float) * C->size);
    }

    while(!stop_criteria){
        // 1. compute the distance matrix
        if (changes) {
            #pragma omp parallel for shared(distance_matrix, C)
            int i;
            for (i = 0; i < C->size; i++) {
                int j;
                for (j = 0; j < C->size; j++) {
                    // TODO: implement distance function
                    distance_matrix[i][j] = distance(C->compressedSets[i], C->compressedSets[j]);
                    distance_matrix[j][i] = distance_matrix[i][j];
                }
            }
        }

        // 3. find the minimum distance in the distance matrix
        float min_distance = FLT_MAX;
        int min_i, min_j;

        #pragma omp parallel for shared(distance_matrix, min_distance, min_i, min_j)
        int i;
        for (i = 0; i < C->size; i++) {
            int j;
            for (j = 0; j < C->size; j++) {
                # pragma omp critical
                if (distance_matrix[i][j] < min_distance) {
                    min_distance = distance_matrix[i][j];
                    min_i = i;
                    min_j = j;
                }
            }
        }

        if (min_distance == FLT_MAX) {
            stop_criteria = true;
            break;
        }

        // 4. merge the two compressed set with the minimum distance
        CompressedSet * newcset = merge_cset(compressedSets[min_i], C->compressedSets[min_j]);

        // 5. check the tightness of the new cluster
        if (tightness(new_cset> T) {
            // 6. add the new cluster to the compressed sets
            add_cset_to_compressed_sets(C, new_cset);
            // 7. remove the two clusters from the compressed sets
            remove_cluster_from_compressed_sets(C, C->compressedSets[min_i]);
            remove_cluster_from_compressed_sets(C, C->compressedSets[min_j]);
            changes = true;
        } else {
            // 8. set the distance between the two clusters to infinity
            distance_matrix[min_i][min_j] = FLT_MAX;
            distance_matrix[min_j][min_i] = FLT_MAX;
            changes = false;

            free(newcset);
        }

    }
}

 load_data(FILE *inputFile, Point *data_buffer, int offset, int rank, int round) {
    // TODO: implement the function
    // load DATA_BUFFER_SIZE points from the input file
    if (stream_cursor == NULL){
        printf("Error: invalid stream cursor\n");
        return false;
    }
}

void add_cluster_to_compressed_sets(CompressedSet *compressedSets, Cluster c) {
    // TODO: copy the function from the serial version trying to multithread it
}

void add_cset_to_compressed_sets(CompressedSet *compressedSets, CompressedSet c) {
    // TODO: copy the function from the serial version trying to multithread it
}

bool primary_compression_criteria(Cluster *clusters, Point p) {
    // function copied from the serial version
    /*
    * Description:
    * Use the Mahalanobis distance to determine whether or not a point 
    * can be added directly to a cluster.
    *
    * Algorithm:
    *   1. check if point is close enough to some centroid, by using mahalanobis radius as a confidence measure
    *   2. if it is, add it to the cluster
    *   3. if not check the secondary compression criteria
    *
    * Parameters:
    *   - clusters: array of clusters
    *   - p: point to check
    *
    * Returns:
    *   - true if point is added to cluster
    *   - false if point is not added to cluster
    * IMPORTANT: the two primary compression criteria are separate and mutually exclusive!
    */

    int i = 0, min_cluster;
    double current_distance, min_distance = DBL_MAX;

    if(DEBUG) printf("      Checking mahalanobis distance for point.\n");
    # pragma omp parallel for shared(min_distance, min_cluster, clusters, p)
    for (i = 0; i < K; i++){
        current_distance = mahalanobis_distance(clusters[i], p);
        // if(DEBUG) printf("      Current distance: %lf.\n", current_distance);
        # pragma omp critical
        if (current_distance < min_distance) {
            min_distance = current_distance;
            min_cluster = i;
        }
    }

    if (min_distance < T){
        if(DEBUG) printf("      Minimal distance is %lf and is under threshold, updating cluster %d with point.\n", min_distance, min_cluster);
        //add point to cluster whose distance from the centroid is minimal, if distance != 0. (the point is not the starting centroid)
        if(min_distance != 0.) update_cluster(&clusters[min_cluster], p); //ALERT: Function implemented in the file "bfr_structure.c", since has costant complexity it remains unvariated
        if(DEBUG) printf("      Cluster %d updated.\n", min_cluster);
        return true;
    }
    
    if(DEBUG) printf("      No cluster fits the criteria, minimal distance is %lf.\n", min_distance);

    return false;
}

void secondary_compression_criteria(Cluster *clusters, RetainedSet *retainedSet, CompressedSet *compressedSets) {
    // TODO: implement the function
    /*
    * Description:
    * find a way to aggregate the retained set into miniclusters,
    * and miniclusters with more than one element can be summarized in the compressed set.
    * Outliers are kept in the retained set.
    *
    * Algorithm:
    *   1. cluster retained set R with classical K-Means, creating k2 clusters
    *   2. clusters that have a tightness measure above a certain threshold are added to the CompressedSet C
    *   3. outliers are kept in the RetainedSet R
    *   4. try aggregating compressed sets using statistics and hierchieal clustering
    *
    * Parameters:
    *   - R: retained set
    *   - clusters: array of clusters
    *   - C: compressed sets
    *
    * Returns:
    *   - void
    */

    // 1. cluster retained set R with classical K-Means, creating k2 clusters, also keep in R the outlayers
    // TODO: implement the function K-Means with open mp
    int k2;
    Cluster *k2_clusters = cluster_retained_set_thrs(retainedSet, k2); //ALERT: not implemented yet

    // 2. clusters that have a tightness measure above a certain threshold are added to the CompressedSet C 
    // using open mp
    #pragma omp parallel for shared(k2_clusters, C)
    for (int i = 0; i < k2; i++){
       add_cluster_to_compressed_sets(C, k2_clusters[i]);
    }

    free(k2_clusters);

    // 4. try aggregating compressed sets using statistics and hierchieal clustering
    hierchieal_clustering_thr(C); //ALERT: not implemented yet
}

bool read_point(Point * data_buffer, Point * p, long int size_of_data_buffer, int * offset){
    //ALERT: FUNCTION COPIED FROM SERIAL VERSION
    
    /*
    * Read a point from data buffer
    *
    * Algorithm:
    *   1. read M coordinates from data buffer
    *   2. if data buffer end is reached return false
    *   3. else return true
    *
    * Parameters:
    *   - data_buffer: buffer to read data from
    *   - p: point to read
    *
    * Returns:
    *   - true if point is read successfully
    *   - false if EOF is reached
    */

    if(DEBUG) printf("  Reading point.\n");
    if (*offset >= size_of_data_buffer){
        return false;
    }

    //read M coordinates from data buffer
    //TODO: check if this is the correct way to read from buffer
    // ---> assuming that the interpretation of the buffer as a static array of point is valid, this should be fixed
    int i = 0;
    for (i = 0; i < M; i++){
        p->coords[i] = (data_buffer[*offset]).coords[i];
    }
    (*offset)++;
    return true;
}

void StreamPoints(Cluster *clusters, CompressedSet *compressedSets, RetainedSet *retainedSet, Point *data_buffer, int size) {
    // TODO: implement the function
    // perform primary compression criteria and secondary compression criteria

    // primary compression criteria
    // for each point in the data buffer
    // TODO: revise this part of the code
    Point p;
    int offset;

    #pragma omp parallel for private(p) shared(data_buffer, clusters, R) reduction(+:offset)
    for (offset = 0; offset < size_of_data_buffer; offset++) {
        if(read_point(data_buffer, &p, size_of_data_buffer, &offset)==0) {
            if (!primary_compression_criteria(clusters, p)) {
                // this function is the same of the serial version, can be found in "bfr_structure.h"
                add_point_to_retained_set(R, p);
            }
        }else {
            break;
        }
    }


    // secondary compression criteria
    secondary_compression_criteria(clusters, R, C);
}

void UpdateCentroids(Cluster *clusters) {
    // TODO: revise the function
    // update the centroids of the clusters
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < DIMENSION; j++) {
            clusters[i].centroid[j] = clusters[i].sum[j] / clusters[i].size;
        }
    }
}


void PlotResults(Cluster *clusters, RetainedSet *retainedSet, CompressedSet *compressedSets) {
    // TODO: implement the function
    // print the results in a file
}
    


int main(int argc, char** argv) {
    int rank, size;

    //set the number of threads
    omp_set_num_threads(NUMBER_OF_THREADS);

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    // catch exceptions
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
        printf("Error: MPI_Comm_rank\n");
        exit(1);
    }

    if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) {
        printf("Error: MPI_Comm_size\n");
        exit(1);
    }

    if (DEBUG) {
        printf("Hello from %d of %d\n", rank, size);
    }

    //initialize the clusters the retained set and the compressed sets
    Cluster *clusters;
    RetainedSet *retainedSet;
    CompressedSet *compressedSets;

    //initialize the data streams from the input files
    FILE *inputFile;
    if (rank == MASTER) {
        // TODO: revise the file opening
        inputFile = fopen("input.txt", "r");
        if (inputFile == NULL) {
            printf("Error opening input file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Broadcast the input file to all processes
        MPI_Bcast(&inputFile, 1, MPI_FILE, MASTER, MPI_COMM_WORLD);
    }else{
        MPI_Recv(&inputFile, 1, MPI_FILE, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    Point * data_buffer; 

    int offset = rank * DATA_BUFFER_SIZE;
    int round = 0;

    //use derived datatype to send the clusters to all processes
    MPI_Datatype arrayclusterType; // array of cluster

    // TODO: implement the derived datatype for the cluster

    // non parallel part of the algorithm
    // init of the clusters made by the master
    if (rank == MASTER) {
        // load the data buffer from the input file with DATA_BUFFER_SIZE points
        bool flag_loaded_data = load_data(inputFile, data_buffer, offset, rank, round);


        // increment the round
        round = round + 1;
        // init the clusters with the centroids
        // can be done in multithread
        // as reference to implement see take_k_centorids() in serial version
        clusters = initClustersWithCentroids(data_buffer, DATA_BUFFER_SIZE);

        
        for (int i = 1; i < size; i++) {
            MPI_Send(clusters, 1, clusterType, i, 0, MPI_COMM_WORLD);
        }
    } else {
        // receive the clusters from the master
        MPI_Recv(clusters, 1, clusterType, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    bool stop_criteria = false;

    // Read the input file simultaneously
    do{
        // update the offset
        offset = rank * DATA_BUFFER_SIZE + round * size * DATA_BUFFER_SIZE;
        // load the data buffer from the input file with DATA_BUFFER_SIZE points
        bool flag_loaded_data = load_data(inputFile, data_buffer, offset, rank, round);
        // increment the round
        round = round + 1;

        // catch exceptions data_buffer is NULL or data_buffer is not full
        // TODO: ALERT THIS IF IS ERRONEOUS
        if (data_buffer == NULL || data_buffer.size() < DATA_BUFFER_SIZE) {
            stop_criteria = true;
        }

        if (!stop_criteria) {
            // perform primary compression criteria and secondary compression criteria
            StreamPoints(clusters, compressedSets, retainedSet, data_buffer, DATA_BUFFER_SIZE);

            // the master process gets the clusters data from the other processes
            // and updates the clusters
            if (rank == MASTER) {
                // the master receive the clusters from the other processes
                // create a temporary clusters to store the clusters from the other processes   
                Cluster *tempClusters;

                tempClusters = (Cluster *)malloc(K * sizeof(Cluster));

                int i = 1;
                for (; i < size; i++) {
                    MPI_Recv(tempClusters, 1, arrayclusterType, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    // update the clusters with the clusters from the other processes
                    // TODO: is important to acknowledge that after the first round the values of the clusters for each dimension should be erased, if not the clusters will be updated readding also the values of the previous rounds but only the master will have the correct values
                    int j = 0;

                    # pragma omp parallel for shared(tempClusters, clusters)
                    for (; j < K; j++) {
                        // TODO: NAIVE IMPLEMENTATION revise
                        clusters[j].size += tempClusters[j].size;
                        int d = 0;
                        for (; d < DIMENSION; d++) {
                            clusters[j].sum[d] += tempClusters[j].sum[d];
                            clusters[j].sum_square[d] += tempClusters[j].sum_square[d];
                        }
                    }
                }
                free(tempClusters);

                // update the clusters
                UpdateCentroids(clusters);
                // send the updated clusters to the other processes
                int i = 1;
                # pragma omp parallel for shared(clusters)
                for (; i < size; i++) {
                    MPI_Send(clusters, 1, arrayclusterType, i, 0, MPI_COMM_WORLD);
                }
            } else {
                // send the clusters to the master
                MPI_Send(clusters, 1, arrayclusterType, MASTER, 0, MPI_COMM_WORLD);

                //create a temporary clusters to store the clusters from the master
                Cluster *tempClusters;

                tempClusters = (Cluster *)malloc(K * sizeof(Cluster));

                // receive the updated clusters from the master
                MPI_Recv(tempClusters,1, arrayclusterType, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // update the clusters centroids with the clusters from the master
                int j = 0;
                for (; j < K; j++) {
                    clusters[j].size += tempClusters[j].size;
                    clusters[j].centroid = tempClusters[j].centroid;
                }

                free(tempClusters);

                // erase the values of the clusters for each dimension in the array sum and sum_square
                int j = 0;
                for (; j < K; j++) {
                    int d = 0;
                    for (; d < DIMENSION; d++) {
                        clusters[j].sum[d] = 0;
                        clusters[j].sum_square[d] = 0;
                    }
                }



            }

            // TODO: we need to decided if the retained set and the compressed sets should be updated in the master or keep in the local processes
            // in the first case we need to send the retained set and the compressed sets to the master and then the master will update the retained set and the compressed sets


            if (rank == MASTER) {
                if(DEBUG) {
                    printf("Round %d\n", round);
                    for (int i = 0; i < K; i++) {
                        printf("Cluster %d: ", i);
                        for (int j = 0; j < DIMENSION; j++) {
                            printf("%f ", clusters[i].centroid[j]);
                        }
                        printf("\n");
                    }
                }
            }

        }
    }while(!stop_criteria);
    
    // Close the input file
    if (rank == MASTER) {
        fclose(inputFile);
    }

    if (rank == MASTER) {
        // print the results
        PlotResults(clusters, retainedSet, compressedSets);

    }

    // free the memory
    free(clusters);
    free(retainedSet);
    free(compressedSets);

    // free the derived datatype
    MPI_Type_free(&arrayclusterType);

    // Finalize the MPI environment
    MPI_Finalize();
    return 0;
}
