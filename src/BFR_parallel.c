// Description: This file contains the parallel implementation of the BFR algorithm.
#include "./BFR_parallel.h"

// TODO: there will some problem here



bool UpdateCentroidsMultithr(Cluster *clusters, int index, int number_of_clusters, int dimension) {
    bool flag_error = false;

    // the function has the same structure of the serial version but in multithreaded version
    int i, j;
    # pragma omp parallel for shared(clusters, index, number_of_clusters, dimension)
    for (i = 0; i < number_of_clusters; i++) {
        for (j = 0; j < dimension; j++) {
            clusters[index].centroid.coords[j] = clusters[i].sum[j] / clusters[i].size;
        }
    }

    return flag_error;
}

Cluster *initClustersWithCentroids(Point *data_buffer, int size, int number_of_clusters, int dimension) {
    // the function has the same structure of the serial version but in multithread version and with some adjustment
    
    // create a set of clusters number_of_clusters
    Cluster *clusters = (Cluster *)malloc(number_of_clusters * sizeof(Cluster));

    // take a random point from the data buffer and make it the centroid of the cluster
    int i, random_index;
    srand(time(NULL));
    random_index = rand() % size;

    //set the index of the cluster
    for(i = 0; i < number_of_clusters; i++) {
        clusters[i].index = i;
    }

    Point centroids[number_of_clusters];

    centroids[0] = data_buffer[random_index];

    int j;
    # pragma omp parallel for shared(clusters, data_buffer, random_index)
    for (j = 0; j < dimension; j++) {
        clusters[0].sum[j] = centroids[0].coords[j];
        clusters[0].sum_squares[j] = centroids[0].coords[j] * centroids[0].coords[j];
    }

    // update the centroids of the clusters
    UpdateCentroidsMultithr(clusters, 0, number_of_clusters, dimension);
    

    // choose the other number_of_clusters-1 centroids seeking the farthest point from the previous centroids
    for (i = 1; i < number_of_clusters; i++) {
        int farthest_point_index = 0;
        float farthest_distance = 0;
        int j;
        // multithread the for loop to find the farthest point from the previous centroids
        // # pragma omp parallel for shared(clusters, data_buffer, farthest_point_index, farthest_distance)
        for (j = 0; j < size; j++) {
            float current_distance = 0;
            int k;
            # pragma omp parallel for shared(clusters, data_buffer, farthest_point_index, farthest_distance, current_distance)
            for (k = 0; k < i; k++) {
                current_distance += distance((Pointer) & centroids[k], (Pointer) & data_buffer[j]);
            }
            // critical section to update the farthest point
            # pragma omp critical
            if (current_distance > farthest_distance) {
                farthest_distance = current_distance;
                farthest_point_index = j;
            }
        }

        // set the farthest point as the centroid of the cluster
        centroids[i] = data_buffer[farthest_point_index];

                // multithread the for loop to update the sum and sum_square of the cluster
        # pragma omp parallel for shared(clusters, data_buffer, farthest_point_index)
        for (j = 0; j < dimension; j++) {
            clusters[i].sum[j] = data_buffer[farthest_point_index].coords[j];
            clusters[i].sum_squares[j] = data_buffer[farthest_point_index].coords[j] * data_buffer[farthest_point_index].coords[j];
        }

        // update the centroids of the clusters in multithreaded version
        UpdateCentroidsMultithr(clusters, i, number_of_clusters, dimension);
    }

    return clusters;
}

CompressedSet *merge_cset(CompressedSet *c1, CompressedSet *c2) {
    // allocate memory for the new compressed set
    CompressedSet * new_cset = (CompressedSet *)malloc(sizeof(CompressedSet));


    // merge the two compressed sets 
    new_cset -> number_of_points = c1 -> number_of_points + c2 -> number_of_points;
    
    int i;
    # pragma omp parallel for shared(new_cset, c1, c2)
    for ( i = 0; i < M; i++){
        new_cset -> sum[i] = c1 -> sum[i] + c2 -> sum[i];
        new_cset -> sum_square[i] = c1 -> sum_square[i] + c2 -> sum_square[i];
    }

    return new_cset;
}

void add_cset_to_compressed_sets(CompressedSets *compressedSets, CompressedSet *c) {
    // TODO: copy the function from the serial version trying to multithread it
}

void hierachical_clustering_thr(CompressedSets * C){
    // TODO: implement this function
    // hierchieal clustering of the compressed sets

    bool stop_criteria = false;
    bool changes = false;

    float **distance_matrix =  malloc (sizeof(float *) * C->number_of_sets);
    int i;
    for (i = 0; i < C->number_of_sets; i++) {
        distance_matrix[i] = malloc (sizeof(float) * C->number_of_sets);
    }

    while(!stop_criteria){
        // 1. compute the distance matrix
        if (changes) {
            int i, j;
            # pragma omp parallel for shared(distance_matrix, C) private(j)
            for (i = 0; i < C->number_of_sets; i++) {
                for (j = 0; j < C->number_of_sets; j++) {
                    // TODO: implement distance function
                    distance_matrix[i][j] = distance((Pointer) & (C->sets[i]), (Pointer) & (C->sets[j]));
                    distance_matrix[j][i] = distance_matrix[i][j];
                }
            }
        }

        // 3. find the minimum distance in the distance matrix
        float min_distance = FLT_MAX;
        int min_i, min_j;

        int i, j;
        # pragma omp parallel for shared(distance_matrix, min_distance, min_i, min_j) private(j)
        for (i = 0; i < C->number_of_sets; i++) {
            for (j = 0; j < C->number_of_sets; j++) {
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
        CompressedSet * new_cset = merge_cset(&(C->sets[min_i]), &(C->sets[min_j]));

        // 5. check the tightness of the new cluster
        if (tightness(new_cset) > T) {
            // 6. add the new cluster to the compressed sets
            add_cset_to_compressed_sets(C, new_cset);
            // 7. remove the two clusters from the compressed sets
            remove_cset_from_compressed_sets(C, C->sets[min_i]);
            remove_cset_from_compressed_sets(C, C->sets[min_j]);
            changes = true;
        } else {
            // 8. set the distance between the two clusters to infinity
            distance_matrix[min_i][min_j] = FLT_MAX;
            distance_matrix[min_j][min_i] = FLT_MAX;
            changes = false;

            free(new_cset);
        }

    }
}

 load_data(FILE *inputFile, Point *data_buffer, int offset, int rank, int round) {
    // TODO: implement the function
    // load DATA_BUFFER_SIZE points from the input file
    if (inputFile == NULL){
        printf("Error: invalid stream cursor\n");
        return false;
    }
}

void add_cluster_to_compressed_sets(CompressedSets *compressedSets, Cluster *c) {
    // TODO: copy the function from the serial version trying to multithread it
}

void add_cset_to_compressed_sets(CompressedSets *compressedSets, CompressedSet * c) {
    compressedSets -> number_of_sets += 1;

    compressedSets -> sets = (CompressedSet *)realloc(compressedSets -> sets, compressedSets -> number_of_sets * sizeof(CompressedSet));

    CompressedSet temp = *c;

    compressedSets -> sets[compressedSets -> number_of_sets - 1] = temp;

    free(c);
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

Cluster *cluster_retained_set_thrs(RetainedSet *R, int *k, int rank, int size){

    if(DEBUG && rank == MASTER) printf("          Initializing standard kmeans data.\n");
    Cluster * miniclusters = init_cluster((*k));

    // TODO: discuss a correct limit for not running standard kmeans
    // as of now, do not run kmeans if the number of points is < k
    // we may want to run kmeans when we have more than k*constant number of points
    if((*R).number_of_points < (*k)){
        if(DEBUG && rank == MASTER) printf("          Retained set has less points (%d) than clusters(%d). Returning empty miniclusters.\n", (*R).number_of_points, (*k));
        return miniclusters;
    }

    kmeans_config config = init_kmeans_config((*k), R, true, rank, size);
    if(DEBUG && rank == MASTER) printf("          Executing standard kmeans.\n");

    if(rank == MASTER){
        kmeans_result result = kmeans(&config);

        if(DEBUG) printf("          Iteration count: %d\n", config.total_iterations);
        if(DEBUG) printf("          Transferring kmeans cluster data to miniclusters.\n");
        
        int i;
        for (i = 0; i < config.num_objs; i++){
            Point *pt = (Point *)(config.objs[i]);

            update_cluster(&miniclusters[config.clusters[i]], *pt);
        }

        // create new correct retained set with only the points left alone in their clusters
        RetainedSet new_R = init_retained_set();
        int * tightness_flag;
        tightness_flag = calloc((*k), sizeof(int));
        for (i = 0; i < config.num_objs; i++){
            Point *pt = (Point *)(config.objs[i]);
            // TODO: use a different measure to determine a minicluster's tightness
            int index = config.clusters[i];

            if (!tightness_evaluation_cluster(miniclusters[index], tightness_flag, index)){
                add_point_to_retained_set(&new_R, *pt);
            }
            else {
                if(DEBUG) printf("Point not added to retained set: %g\t%g\t%d\n", pt->coords[0], pt->coords[1], config.clusters[i]);
            }
        }

        if(DEBUG){
            printf("          Old retained set:\n");
            print_retainedset(*R);
        }

        // TODO: discuss whether this is correct or not
        // free old retained set and replace with new one
        free((*R).points);
        (*R).points = new_R.points;
        (*R).number_of_points = new_R.number_of_points;

        if(DEBUG){
            printf("          New retained set:\n");
            print_retainedset(*R);
        }

        if(DEBUG) printf("          Freeing previously allocated data for standard kmeans.\n");

        //update miniclusters, retain only with tightness_flag = 2
        (*k) = update_miniclusters(&miniclusters, tightness_flag, (*k));

        // free the kmeans' config data
        free(config.objs);
        free(config.centers);
        free(config.clusters);

        // free tightness flag
        free(tightness_flag);

        // update miniclusters' centroids
        update_centroids(&miniclusters, (*k));

    }
    else{
        kmeans(&config);
        free(config.objs);
        free(config.centers);
        free(config.clusters);
    }

    // If MASTER, will be correctly filled. If not, it's an empty array
    return miniclusters;
}

void secondary_compression_criteria(Cluster *clusters, RetainedSet *retainedSet, CompressedSets *compressedSets, int rank, int size) {
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
    int k2 = K;
    Cluster *k2_clusters = cluster_retained_set_thrs(retainedSet, &k2, rank, size);

    // 2. clusters that have a tightness measure above a certain threshold are added to the CompressedSet C 
    // using open mp
    int i;
    # pragma omp parallel for shared(k2_clusters, clusters)
    for (i = 0; i < k2; i++){
       //    add_cluster_to_compressed_sets(clusters, k2_clusters[i]);
    }

    free(k2_clusters);

    // 4. try aggregating compressed sets using statistics and hierachical clustering
    // hierachical_clustering_thr(clusters); //ALERT: not implemented yet
}

void UpdateRetainedSet(RetainedSet *R, RetainedSet *tempRetainedSet){
    int old_number_of_points = R->number_of_points;
    R->number_of_points += tempRetainedSet->number_of_points;
    R->points = realloc(R->points, R->number_of_points * sizeof(Point));
    if (R->points == NULL){
        perror("Error: could not allocate memory\n");
        exit(1);
    }

    int i;
    for (i = 0; i < tempRetainedSet->number_of_points; i++){
        R->points[old_number_of_points - 1 + i] = tempRetainedSet->points[i];
    }
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

void StreamPoints(Cluster *clusters, CompressedSets *compressedSets, RetainedSet *retainedSet, Point *data_buffer, int size) {
    // TODO: implement the function
    // perform primary compression criteria and secondary compression criteria

    // primary compression criteria
    // for each point in the data buffer
    // TODO: revise this part of the code
    Point p;
    int offset;

    # pragma omp parallel for private(p) shared(data_buffer, clusters, retainedSet) //reduction(+:offset)
    for (offset = 0; offset < size; offset++) {
        if(read_point(data_buffer, &p, size, &offset)==0) {
            if (!primary_compression_criteria(clusters, p)) {
                // this function is the same of the serial version, can be found in "bfr_structure.h"
                add_point_to_retained_set(retainedSet, p);
            }
        }else {
            // break;
            continue;
        }
    }
}


void PlotResults(Cluster *clusters, RetainedSet *retainedSet, CompressedSets *compressedSets) {
    // TODO: implement the function
    // print the results in a file
}

void remove_cset_from_compressed_sets(CompressedSet *C, CompressedSet * c) {
    C->number_of_sets -= 1;

    CompressedSets * temp = (CompressedSet *)malloc(sizeof(CompressedSet) * C->number_of_sets);

    int i;
    for (i = 0; i < C->number_of_sets; i++) {
        if (C->sets[i] == c) {
            continue;
        }
        temp[i] = C->sets[i];
    }

    free(C->sets);
    free(C);

    C = temp;
}

bool hierachical_clust_parallel(CompressedSets * C, int rank, int size){
    /*
    * Description:
    * hierchieal clustering of the compressed sets
    *
    * Algorithm:
    *   1. compute the distance matrix
    *   2. find the minimum distance in the distance matrix
    *   3. merge the two compressed set with the minimum distance
    *   4. check the tightness of the new cluster
    *   5. if the new cluster is tight enough add it to the compressed sets
    *   6. remove the two clusters from the compressed sets
    *   7. else set the distance between the two clusters to infinity
    *   8. repeat from step 1 until the stop criteria is reached
    *
    * Parameters:
    *   - C: compressed sets
    *
    * Returns:
    *   - void
    */
    // initialize the distance matrix
    float **distance_matrix =  malloc (sizeof(float *) * C->number_of_sets);
    int i;
    for (i = 0; i < C->number_of_sets; i++) {
        distance_matrix[i] = malloc (sizeof(float) * C->number_of_sets);
    }

    bool stop_criteria = false;
    bool changes = false;
    int iterations = 0;


    while(!stop_criteria){
        int start, number_of_sets;

        // 1. compute the distance matrix
        if (changes) {
            // free the distance matrix
            for (i = 0; i < C->number_of_sets; i++) {
                free(distance_matrix[i]);
            }
            free(distance_matrix);

            distance_matrix =  malloc (sizeof(float *) * C->number_of_sets);

            for (i = 0; i < C->number_of_sets; i++) {
                distance_matrix[i] = malloc (sizeof(float) * C->number_of_sets);
            }

            if (size < C->number_of_sets ){
                // case where at least one process has to compute the distance vector more than once
                number_of_sets = C->number_of_sets / size;

                if (rank == size - 1){
                    if (number_of_sets * size < C->number_of_sets){
                        number_of_sets = number_of_sets + C->number_of_sets % size;
                    }else{
                        number_of_sets = C->number_of_sets - (size - 1) * number_of_sets;
                    }
                }

                start = rank * number_of_sets; 
            }else if (size >= C->number_of_sets){
                // case where each process has to compute the distance vector only once
                if (rank > C->number_of_sets){
                    stop_criteria = true;
                    break;
                }

                number_of_sets = 1;
                start = rank;
            }

            // for each process compute the distance vectors for the sets assigned
            int i, j;
            # pragma omp parallel for shared(distance_matrix, C, start, number_of_sets) private(j)
            for (i = start; i < start + number_of_sets; i++) {
                for (j = i; j < C->number_of_sets; j++) {
                    distance_matrix[i][j] = distance((Pointer) & (C->sets[i]), (Pointer) & (C->sets[j]));
                }
            }

            if (rank == MASTER) {
                // case where the master has to compute the distance vector for the remaining sets
                for (i = 1; i < C->number_of_sets; i++) {
                    // the rank of the process that should have computed the distance vector for the set i: i / number_of_points - 1 
                    // if the rank is the master the distance vector is already computed
                    // example: if the number of sets is 10 and the number of processes is 4, the first 3 processes will compute the distance vector for 3 sets, the last process will compute the distance vector for 1 set
                    // sets i would have been computed by the process i / number_of_points - 1, in the example the set 4 would have been computed by the process 1
                    int expected_rank = floor(i / number_of_sets);
                    MPI_Recv(distance_matrix[i], C->number_of_sets, MPI_FLOAT, expected_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }

                //fill the upper part of the distance matrix
                //ALERT: REVISE THIS PART
                for (i = 0; i < C->number_of_sets; i++) {
                    for (j = 0; j < i; j++) {
                        distance_matrix[i][j] = distance_matrix[j][i];
                    }
                }

                //Broadcast the distance matrix to the other processes
                for(i = 1; i < size; i++) {
                    MPI_Send(distance_matrix[i], number_of_sets, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
                }
            } else {
                for(i = start; i < start + number_of_sets; i++) {
                    int size_of_dist_vector = C->number_of_sets-i;
                    MPI_Send(distance_matrix[start], number_of_sets, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD);
                }

                //Receive the distance matrix from the master
                for (i = 0; i < C->number_of_sets; i++) {
                    MPI_Recv(distance_matrix[i], C->number_of_sets, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }            
        }

        // 3. find the minimum distance in the distance matrix

        if (size < C->number_of_sets ){
            // case where at least one process has to compute the distance vector more than once
            number_of_sets = C->number_of_sets / size;

            if (rank == size - 1){
                if (number_of_sets * size < C->number_of_sets){
                    number_of_sets = number_of_sets + C->number_of_sets % size;
                }else{
                    number_of_sets = C->number_of_sets - (size - 1) * number_of_sets;
                }
            }

            start = rank * number_of_sets;
        }else if (size >= C->number_of_sets){
            // case where each process has to compute the distance vector only once
            if (rank > C->number_of_sets){
                stop_criteria = true;
                break;
            }

            number_of_sets = 1;
            start = rank;
        }

        float min_distance = FLT_MAX;
        int min_i, min_j;

        int i, j;
        # pragma omp parallel for shared(distance_matrix, min_distance, min_i, min_j) private(j)
        for (i = start; i < start + number_of_sets; i++) {
            for (j = i; j < C->number_of_sets; j++) {
                # pragma omp critical
                if (distance_matrix[i][j] < min_distance) {
                    min_distance = distance_matrix[i][j];
                    min_i = i;
                    min_j = j;
                }
            }
        }

        // the master process gets the minimum distance from the other processes
        if (rank == MASTER) {
            float temp_min_distance;
            int temp_min_i, temp_min_j;

            for (i = 1; i < size; i++) {
                MPI_Recv(&temp_min_distance, 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&temp_min_i, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&temp_min_j, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (temp_min_distance < min_distance) {
                    min_distance = temp_min_distance;
                    min_i = temp_min_i;
                    min_j = temp_min_j;
                }
            }
        } else {
            MPI_Send(&min_distance, 1, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD);
            MPI_Send(&min_i, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
            MPI_Send(&min_j, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
        }
        
        if (min_distance >= FLT_MAX - 100. || iterations > UPPER_BOUND_ITERATIONS) {
            stop_criteria = true;
            break;
        }


        // 4. merge the two compressed set with the minimum distance in the master
        if (rank == MASTER){
            CompressedSet * new_cset = merge_cset(&(C->sets[min_i]), &(C->sets[min_j]));

            // 5. check the tightness of the new cluster
            if (tightness(new_cset) > T) {
                // 6. add the new cluster to the compressed sets
                CompressedSet * temp = (CompressedSet *)malloc(sizeof(CompressedSet));
                *temp = *new_cset;

                add_cset_to_compressed_sets(C, temp);
                // 7. remove the two clusters from the compressed sets
                remove_cset_from_compressed_sets(C, &C->sets[min_i]);
                remove_cset_from_compressed_sets(C, &C->sets[min_j]);
                changes = true;

                MPI_Bcast(new_cset->number_of_points, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
                MPI_Bcast(new_cset->sum, M, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
                MPI_Bcast(new_cset->sum_square, M, MPI_FLOAT, MASTER, MPI_COMM_WORLD);

                free(new_cset);
            } else {
                // 8. set the distance between the two clusters to infinity
                distance_matrix[min_i][min_j] = FLT_MAX;
                distance_matrix[min_j][min_i] = FLT_MAX;
                changes = false;

                free(new_cset);
            }

            //Broadcast the changes to the other processes
            MPI_Bcast(&changes, 1, MPI_C_BOOL, MASTER, MPI_COMM_WORLD);
            MPI_Bcast(&min_i, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
            MPI_Bcast(&min_j, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        }else{
            MPI_Recv(&changes, 1, MPI_C_BOOL, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&min_i, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&min_j, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

             if (changes) {
                CompressedSet * new_cset = (CompressedSet *)malloc(sizeof(CompressedSet));
                int number_of_points;
                float sum[M], sum_square[M];

                MPI_Recv(&number_of_points, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(sum, M, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(sum_square, M, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                new_cset -> number_of_points = number_of_points;
                int i;
                # pragma omp parallel for shared(new_cset, sum, sum_square)
                for (i = 0; i < M; i++){
                    new_cset -> sum[i] = sum[i];
                    new_cset -> sum_square[i] = sum_square[i];
                }

                add_cset_to_compressed_sets(C, new_cset);
                remove_cset_from_compressed_sets(C, &C->sets[min_i]);
                remove_cset_from_compressed_sets(C, &C->sets[min_j]);

            }else{
                distance_matrix[min_i][min_j] = FLT_MAX;
                distance_matrix[min_j][min_i] = FLT_MAX;
            }
        }

        iterations++;
    }
    

    // synchronize the processes
    MPI_Barrier(MPI_COMM_WORLD);

    // free the distance matrix
    for (i = 0; i < C->number_of_sets; i++) {
        free(distance_matrix[i]);
    }
    free(distance_matrix);

    return true;
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
    RetainedSet retainedSet_normal = (init_retained_set());
    CompressedSets compressedSets_normal = (init_compressed_sets());
    RetainedSet *retainedSet = &retainedSet_normal;
    CompressedSets *compressedSets = &compressedSets_normal;

    //init retained set and compressed sets
    // retainedSet = (RetainedSet *)malloc(sizeof(RetainedSet));
    // retainedSet -> number_of_points = 0;
    // compressedSets = (CompressedSet *)malloc(sizeof(CompressedSet));
    // compressedSets -> number_of_sets = 0;

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
        MPI_Bcast(&inputFile, 1, MPI_File, MASTER, MPI_COMM_WORLD);
    }else{
        MPI_Recv(&inputFile, 1, MPI_File, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    Point * data_buffer; 

    int offset = rank * DATA_BUFFER_SIZE;
    int round = 0;

    //use derived datatype to send the clusters to all processes
    MPI_Datatype arrayclusterType; // array of cluster
    MPI_Datatype MPI_RETAINED_SET;

    // TODO: implement the derived datatype for the cluster and for the Retained Set

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
        clusters = initClustersWithCentroids(data_buffer, DATA_BUFFER_SIZE, K, DIMENSION);

        int i;        
        for (i = 1; i < size; i++) {
            //TODO: adjust this send and next receive
            MPI_Send(clusters, 1, arrayclusterType, i, 0, MPI_COMM_WORLD);
        }
    } else {
        // receive the clusters from the master
        MPI_Recv(clusters, 1, arrayclusterType, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
            // perform primary compression criteria
            StreamPoints(clusters, retainedSet, data_buffer, DATA_BUFFER_SIZE);

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
                    // TODO: is important to acknowledge that after the first round the values of the clusters for each M should be erased, if not the clusters will be updated readding also the values of the previous rounds but only the master will have the correct values
                    int j, d;

                    # pragma omp parallel for shared(tempClusters, clusters) private(d)
                    for (j = 0; j < K; j++) {
                        // TODO: NAIVE IMPLEMENTATION revise
                        clusters[j].size += tempClusters[j].size;
                        for (d = 0; d < M; d++) {
                            clusters[j].sum[d] += tempClusters[j].sum[d];
                            clusters[j].sum_squares[d] += tempClusters[j].sum_squares[d];
                        }
                    }
                }
                free(tempClusters);
                for (i = 0; i < K; i++) {
                    // update the clusters
                    UpdateCentroidsMultithr(clusters, 0, K, DIMENSION);
                }
                // send the updated clusters to the other processes
                
                # pragma omp parallel for shared(clusters)
                for (i = 1; i < size; i++) {
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

                // erase the values of the clusters for each M in the array sum and sum_square
                j = 0;
                for (; j < K; j++) {
                    int d = 0;
                    for (; d < M; d++) {
                        clusters[j].sum[d] = 0;
                        clusters[j].sum_squares[d] = 0;
                    }
                }
            }

            // TODO: we need to decided if the retained set and the compressed sets should be updated in the master or keep in the local processes
            // in the first case we need to send the retained set and the compressed sets to the master and then the master will update the retained set and the compressed sets

            // Reducing the retained set, by gathering the retained set from the other processes and updating the retained set
            if (rank == MASTER) {
                // reduce the retained set
                RetainedSet *tempRetainedSet;
                tempRetainedSet = (RetainedSet *)malloc(sizeof(RetainedSet));

                int i = 1;
                for (; i < size; i++) {
                    MPI_Recv(tempRetainedSet, 1, MPI_RETAINED_SET, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    // update the retained set with the retained set from the other processes
                    UpdateRetainedSet(retainedSet, tempRetainedSet);
                }

                free(tempRetainedSet);
            } else {
                // send the retained set to the master
                MPI_Send(retainedSet, 1, MPI_RETAINED_SET, MASTER, 0, MPI_COMM_WORLD);

                free(retainedSet->points);
                free(retainedSet);

                retainedSet = (RetainedSet *)malloc(sizeof(RetainedSet));
                retainedSet->number_of_points = 0;
            }

            // TODO: perform the secondary compression criteria in parallel version
            secondary_compression_criteria(clusters, retainedSet, compressedSets, rank, size);

            // TODO: insert the kmeans clustering in the parallel version, implemented in "kmeans.c"

            // TODO: filtering of cluster in the master by tightness criterion

            // TODO: merge the clusters in the compressed sets in the master

            // TODO: perform the hierachical clustering in parallel version

            hierachical_clust_parallel(compressedSets);


            // TODO: synchronize the processes

            if (rank == MASTER) {
                if(DEBUG) {
                    printf("Round %d\n", round);
                    int i, j;
                    for (i = 0; i < K; i++) {
                        printf("Cluster %d: ", i);
                        for (j = 0; j < M; j++) {
                            printf("%f ", clusters[i].centroid.coords[j]);
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
        // TODO: try merge compressed sets with clusters in the master, using tightness criterion evaluate if the clusters can be merged
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
