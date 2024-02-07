// Description: This file contains the parallel implementation of the BFR algorithm.
#include <mpi.h>
#include <omp.h>
#include "../lib/kmeans_wrapper/kmeans_wrapper.h"

// TODO: there will some problem here

#define MASTER 0
#define DEBUG 1
# define K 10
# define DIMENSION 2
# define NUMBER_OF_THREADS 4

cluster *initClustersWithCentroids(FILE *inputFile, Point *data_buffer, int size) {
    // TODO: implement the function, THIS IS JUST A NAIVE IMPLEMENTATION
    
    // take the first k points as centroids
    Cluster *clusters = (Cluster *)malloc(K * sizeof(Cluster));
    for (int i = 0; i < K; i++) {
        clusters[i].centroid = data_buffer[i];
        clusters[i].size = 0;
    }
    return clusters;
}

void load_data(FILE *inputFile, Point *data_buffer, int offset, int rank, int round) {
    // TODO: implement the function
    // load DATA_BUFFER_SIZE points from the input file
    // the offset is the position in the file to start reading
    // rank is the rank of the process
    // round is the round of the algorithm
}

void add_point_to_retained_set(RetainedSet *retainedSet, Point p) {
    // TODO: copy the function from the serial version
}

void primary_compression_criteria(Cluster *clusters, Point p) {
    // TODO: copy the function from the serial version
}

void secondary_compression_criteria(Cluster *clusters, RetainedSet *retainedSet, CompressedSet *compressedSets) {
    // TODO: implement the function
    // perform secondary compression criteria
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
        inputFile = fopen("input.txt", "r");
        if (inputFile == NULL) {
            printf("Error opening input file\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Broadcast the input file to all processes
    MPI_Bcast(&inputFile, 1, MPI_FILE, MASTER, MPI_COMM_WORLD);

    Point data_buffer[DATA_BUFFER_SIZE];

    int offset = rank * DATA_BUFFER_SIZE;
    int round = 0;

    //use derived datatype to send the clusters to all processes
    MPI_Datatype clusterType;

    // TODO: implement the derived datatype for the cluster

    // non parallel part of the algorithm
    // init of the clusters made by the master
    if (rank == MASTER) {
        // load the data buffer from the input file with DATA_BUFFER_SIZE points
        data_buffer = load_data(inputFile, offset, rank, round);
        // increment the round
        round = round + 1;
        // init the clusters with the centroids
        // can be done in multithread
        // as reference to implement see take_k_centorids() in serial version
        clusters = initClustersWithCentroids(inputFile, data_buffer, DATA_BUFFER_SIZE);

        
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
        data_buffer = load_data(inputFile, offset, rank, round);
        // increment the round
        round = round + 1;

        // catch exceptions data_buffer is NULL or data_buffer is not full
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
                int i = 1;
                for (; i < size; i++) {
                    MPI_Recv(tempClusters, 1, clusterType, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    // update the clusters with the clusters from the other processes
                    // TODO: is important to acknowledge that after the first round the values of the clusters for each dimension should be erased, if not the clusters will be updated readding also the values of the previous rounds
                    int j = 0;
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
                // update the clusters
                UpdateCentroids(clusters);
                // send the updated clusters to the other processes
                int i = 1;
                for (; i < size; i++) {
                    MPI_Send(clusters, 1, clusterType, i, 0, MPI_COMM_WORLD);
                }
            } else {
                // send the compressed sets to the master
                MPI_Send(clusters, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
                // receive the updated clusters from the master
                MPI_Recv(clusters, 1, clusterType, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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

        // free the memory
        free(clusters);
        free(retainedSet);
        free(compressedSets);

        // free the derived datatype
        MPI_Type_free(&clusterType);
    }

    // Finalize the MPI environment
    MPI_Finalize();
    return 0;
}
