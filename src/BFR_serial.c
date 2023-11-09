#include "../lib/kmeans_wrapper/kmeans_wrapper.h"

void take_k_centroids(Cluster *clusters, Point * data_buffer, long size_of_data_buffer){
    /*
    * Take K random centroids from input space
    *
    * Algorithm:
    *   1. take K random points from input space
    *   2. set these points as centroids
    *
    * Parameters:
    *   - clusters: array of clusters
    *   - input_file: file to read points from
    *
    * Returns:
    *   - void
    */
    int i = 0;
    srand(time(NULL));
    for (i = 0; i < K; i++){
        // centroid are random points in N^M space
        Point p;
        int j = 0;
        int r = lround(size_of_data_buffer * (1.0 * rand() / RAND_MAX));
        for (j = 0; j < M; j++){
            /* Pointers to the randomly picked point */
            
            p.coords[j] = data_buffer[r].coords[j];
            clusters[i].sum[j] += p.coords[j];
            clusters[i].sum_squares[j] += pow(p.coords[j], 2);
        }
        p.cluster = i;
        clusters[i].centroid = p;
        clusters[i].size = 1;
    }
}

bool read_point(Point * data_buffer, Point * p, long int size_of_data_buffer, int * offset){
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

void merge_points(Point * p1, Point * p2, CompressedSets * C, RetainedSet * R){
    /*
    * Merge two points and create a new compressed set
    *
    * Algorithm:
    *   1. create new compressed set
    *   2. add points to new compressed set
    *   3. add new compressed set to compressed sets
    *   4. remove point from retained set
    *
    * Parameters:
    *   - p1: first point
    *   - p2: second point
    *   - C: compressed sets
    *   - R: retained set
    *
    * Returns:
    *   - void
    */

    //create new compressed set
    CompressedSet set;
    set.number_of_points = 2;
    int j = 0;
    for (j = 0; j < M; j++){
        set.sum[j] = p1->coords[j] + p2->coords[j];
        set.sum_square[j] = pow(p1->coords[j], 2) + pow(p2->coords[j], 2);
    }

    //add new compressed set to compressed sets
    (*C).number_of_sets += 1;
    (*C).sets = realloc(C->sets, C->number_of_sets * sizeof(CompressedSet));
    if (C->sets == NULL){
        printf("Error: could not allocate memory\n");
        exit(1);
    }
    (*C).sets[C->number_of_sets - 1] = set;

    //remove point from retained set
    int size_of_R = R->number_of_points;
    int index = 0;
    for (j = 0; j < size_of_R; j++){
        if (R->points[j].coords[0] == p1->coords[0] && R->points[j].coords[1] == p1->coords[1]){
            index = j;
        }
    }

    for (j = index; j < size_of_R - 1; j++){
        R->points[j] = R->points[j + 1];
    }

    R->number_of_points -= 1;
    R->points = realloc(R->points, R->number_of_points * sizeof(Point));
}

Point perturbate_centroid(Cluster cluster, Point p, bool is_first){
    /*
    * Perturbate centroid
    *
    * Algorithm:
    *   1. calculate confidence interval
    *   2. perturbate centroid
    *
    * Parameters:
    *   - cluster: cluster to perturbate centroid
    *   - p: point to check
    *   - is_first: flag to check if it is the first centroid
    *
    * Returns:
    *   - perturbated centroid
    */

    //TODO: the calculation of confidence interval is not decided, this is just a proposition
    double confidence_interval = 0;
    int i = 0;
    for (i = 0; i < M; i++){
        double mean = (double) cluster.sum[i] / (double) cluster.size;
        double variance = (double) cluster.sum_squares[i] / (double) cluster.size - pow(mean, 2);
        double standard_deviation = sqrt(variance);
        confidence_interval += pow(standard_deviation, 2);
    }
    confidence_interval = sqrt(confidence_interval);

    //TODO: the perturbation of centroid is not decided, this is just a prototype, it has to changed
    Point pertubated_centroid;
    if (is_first){
        pertubated_centroid = cluster.centroid;
        int j = 0;
        for (j = 0; j < M; j++){
            pertubated_centroid.coords[j] += confidence_interval;
        }
    }else{
        pertubated_centroid = cluster.centroid;
        int j = 0;
        for (j = 0; j < M; j++){
            pertubated_centroid.coords[j] -= confidence_interval;
        }
    }

    return pertubated_centroid;
}

bool primary_compression_criteria(Cluster * clusters, Point p){
    /*
    * Idea:
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
    for (i = 0; i < K; i++){
        current_distance = mahalanobis_distance(clusters[i], p);
        // if(DEBUG) printf("      Current distance: %lf.\n", current_distance);
        if (current_distance < min_distance) {
            min_distance = current_distance;
            min_cluster = i;
        }
    }

    if (min_distance < T){
        if(DEBUG) printf("      Minimal distance is %lf and is under threshold, updating cluster %d with point.\n", min_distance, min_cluster);
        //add point to cluster whose distance from the centroid is minimal
        update_cluster(&clusters[min_cluster], p);
        if(DEBUG) printf("      Cluster %d updated.\n", min_cluster);
        return true;
    }
    
    if(DEBUG) printf("      No cluster fits the criteria, minimal distance is %lf.\n", min_distance);

    return false;

}

bool second_primary_compression_criteria(Cluster * clusters, Point p){
    /*
    * Idea:
    * Takes the two closest centroids 
    * if after pertubation of centroids in confidence interval in purpose to set first centroid most far away from the point
    * and the second centroid most close to the point
    * if the point is still close to the first centroid, add it to the cluster
    *
    * Algorithm:
    *   1. take the two closest centroids
    *   2. if after pertubation of centroids in confidence interval in purpose to set first centroid most far away from the point
    *   3. and the second centroid most close to the point
    *   4. if the point is still close to the first centroid, add it to the cluster
    *   5. if not check the secondary compression criteria
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

    //take the two closest centroids
    //TODO: the confidence interval is not correct
    double min_distance = distance((Pointer)&(clusters[0].centroid), (Pointer)&p);
    int index_of_min = 0;
    double second_min_distance = distance((Pointer)&(clusters[1].centroid), (Pointer)&p);
    int index_of_second_min = 1;
    int i;

    for (i = 1; i < K; i++){
        double dist = distance((Pointer)&(clusters[i].centroid), (Pointer)&p);
        if (dist < min_distance){
            second_min_distance = min_distance;
            index_of_second_min = index_of_min;
            min_distance = dist;
            index_of_min = i;
        }else if (dist < second_min_distance){
            second_min_distance = dist;
            index_of_second_min = i;
        }
    }

    Point pertubated_centroid_1 = perturbate_centroid(clusters[index_of_min], p, true); 
    Point pertubated_centroid_2 = perturbate_centroid(clusters[index_of_second_min], p, false);

    if (distance((Pointer)&pertubated_centroid_1, (Pointer)&p) < distance((Pointer)&pertubated_centroid_2, (Pointer)&p)){
        //add point to cluster
        if(DEBUG) printf(" Point respect the second primary compression criteria, adding to cluster.\n");
        update_cluster(&clusters[index_of_min], p);
        if(DEBUG) printf(" Cluster %d updated.\n", index_of_min);
        return true;
    }

    return false;
}

void secondary_compression_criteria(RetainedSet * R, Cluster * clusters, CompressedSets * C){
    /*
    * Idea:
    * find a way to aggregate the retained set into miniclusters,
    * and miniclusters with more than one element can be summarized in the compressed set.
    * Outliers are kept in the retained set.
    *
    * Algorithm:
    *   1. cluster retained set R with classical K-Means, creating k2 clusters
    *   2. clusters with more than one point in it are added to CompressedSets C
    *   3. outliers are kept in the RetainedSet R
    *   4. try aggregating compressed sets using statistics
    *
    * Parameters:
    *   - R: retained set
    *   - clusters: array of clusters
    *   - C: compressed sets
    *
    * Returns:
    *   - void
    */

    // aggregate retained set into miniclusters
    if(DEBUG) printf("      Creating miniclusters from Retained Set.\n");
    // TODO: decide the number of clusters. for now, set as K but by specs should be k2>K
    int number_of_miniclusters = K;
    Cluster *miniclusters = cluster_retained_set(R, &number_of_miniclusters);

    if(DEBUG) printf("      Merging Compressed Sets and miniclusters and aggregating them, if possible.\n");
    merge_compressedsets_and_miniclusters(C, miniclusters, number_of_miniclusters);


    if(DEBUG) printf("      Freeing miniclusters.\n");
    free(miniclusters);

    // if (C->number_of_sets > 0){
    //     int size_of_C = C->number_of_sets;
    //     int j = 0;
    //     for (j = 0; j < size_of_C; j++){
    //         // calculate center of set
    //         Point center;
    //         int k = 0;
    //         for (k = 0; k < M; k++){
    //             center.coords[k] = C->sets[j].sum[k] / C->sets[j].number_of_points;
    //         }
    //         if (distance(p, center) < T){
    //             // add point to set
    //             C->sets[j].number_of_points += 1;
    //             int k = 0;
    //             for (k = 0; k < M; k++){
    //                 C->sets[j].sum[k] += p.coords[k];
    //                 C->sets[j].sum_square[k] += pow(p.coords[k], 2);
    //             }
    //             break;
    //         }
    //     }
    // }
    // else{
    //     // check if point is close to any point in retained set
    //     //  in this case merge the two points and create a new compressed set
    //     bool is_close = false;
    //     int size_of_R = R->number_of_points;
    //     int j = 0;
    //     for (j = 0; j < size_of_R; j++){
    //         if (distance(p, R->points[j]) < T){
    //             // merge points
    //             merge_points(&p, &R->points[j], C, R);
    //             is_close = true;
    //             break;
    //         }
    //     }

    //     if (!is_close){
    //         // add point to retained set
    //         (*R).number_of_points += 1;
    //         (*R).points = realloc(R->points, R->number_of_points * sizeof(Point));
    //         if (R->points == NULL){
    //             printf("Error: could not allocate memory\n");
    //             exit(1);
    //         }
    //         (*R).points[R->number_of_points - 1] = p;
    //     }
    // }
}

void stream_data(RetainedSet * R, Cluster * clusters, CompressedSets * C, Point * data_buffer, long int size_of_data_buffer){
    /*
    * Algorithm:
    *   1. read point from buffer
    *   2. find if point is close enough to some centroid through primary compression criteria
    *   3. if it is, add it to the cluster
    *   4. if not, add it to retained set
    *   5. apply secondary compression criteria to retained set
    * 
    * Parameters:
    *   - R: retained set
    *   - clusters: array of clusters
    *   - C: compressed sets
    *   - data_buffer: buffer to read data from
    *   - size_of_data_buffer: size of data buffer
    *
    * Returns:
    *   - void
    */

    if(DEBUG) printf("  Streaming data.\n");

    Point p;        // current point being read
    int offset = 0; // counter of the # of points read

    if(DEBUG) printf("  Reading single points from buffer.\n");
    while (read_point(data_buffer, &p, size_of_data_buffer, &offset)){
        // apply primary compression criteria
        if(DEBUG) printf("  Checking primary compression criteria.\n");

        if (!primary_compression_criteria(clusters, p)){
            if(DEBUG) printf("  Primary compression criteria not applicable, adding to Retained Set.\n");
            // if not, add point to retained set
            add_point_to_retained_set(R, p);
        }
    }

    if(DEBUG) printf("  Checking secondary compression criteria.\n");
    // apply secondary compression criteria over retained set and compressed sets
    secondary_compression_criteria(R, clusters, C);
}

bool load_data_buffer(data_streamer stream_cursor, Point * data_buffer, int max_size_of_buffer, long * size_of_data_buffer){
    if (stream_cursor == NULL){
        printf("Error: invalid stream cursor\n");
        exit(1);
    }

    // if ((void*)data_buffer == NULL){
    //     data_buffer = (Point *)malloc(max_size_of_buffer * sizeof(Point));
    //     if ((void*)data_buffer == NULL){
    //         printf("Error: could not allocate memory\n");
    //         return false;
    //     }
    // }
    
    // fills data_buffer with MAX_SIZE_OF_BUFFER points or until EOF is reached
    int i, j, counter = 0, return_val;
    float x;
    
    if(DEBUG) printf("  Reading points.\n");
    for (i = 0; i < MAX_SIZE_OF_BUFFER; i++){
        Point p;
        for (j = 0; j < M; j++){
            return_val = fscanf(stream_cursor, "%f", &x);

            if (return_val == EOF && i==0 && j==0) {
                if(DEBUG) printf("  Buffer is empty, end algorithm.\n");
                // first iteration of the loops, buffer is empty, end algorithm
                return true;
            }
            else if (return_val == EOF) {
                if(DEBUG) printf("  Buffer is not empty, end algorithm at the next iteration.\n");
                // buffer is not empty - do another iteration, end at the next one.
                return false;
            }

            p.coords[j] = (double) x;
            data_buffer[i] = p;
        }
        counter++;
    }

    *size_of_data_buffer = counter;
    if(DEBUG) printf("  Data buffer size is %d.\n", counter);
    // *size_of_data_buffer = fread(*data_buffer, 1, max_size_of_buffer, stream_cursor);
    // if (*size_of_data_buffer == 0){
    //     return true;
    // }

    //TODO: check if the size of data respect the point format: n floats separated by space
    // ---> now data should be read regardless of format

    return false;
}

int main(int argc, char ** argv){
    if (argc != 2){
        printf("Usage: %s <input_file>\n", argv[0]);
        return 1;
    }

    if(DEBUG) printf("Starting BFR.\n");

    data_streamer stream_cursor = data_streamer_Init(argv[1], "r");

    Cluster * clusters = init_cluster(K); // allocate the clusters array dynamically

    RetainedSet R = init_retained_set();

    CompressedSets C = init_compressed_sets();

    Point data_buffer[MAX_SIZE_OF_BUFFER];

    if(DEBUG) printf("Allocation phase complete. Starting iterations.\n");

    //start reading block points from input file
    bool stop_criteria = false, first_round = true;
    do{ 
        long size_of_data_buffer = 0;
        if(DEBUG) printf("Loading buffer.\n");
        stop_criteria = load_data_buffer(stream_cursor, data_buffer, MAX_SIZE_OF_BUFFER, &size_of_data_buffer);

        if(DEBUG) printf("Buffer loaded, stop criteria is %d\n", stop_criteria);
        if(stop_criteria == false){

            if(first_round) take_k_centroids(clusters, data_buffer, size_of_data_buffer);

            // printf("data buffer: %s\n", data_buffer);
            //stream data
            //TODO: proposition to change the name of this function, in order to make it more clear
            if(DEBUG) printf("Streaming data.\n");
            stream_data(&R, clusters, &C, data_buffer, size_of_data_buffer);

            if(DEBUG) printf("Updating centroids.\n");
            //update centroids
            update_centroids(&clusters, K);

            if(DEBUG) printf("Printing.\n");
            //print clusters
            print_clusters(clusters);

            //print compressed sets
            print_compressedsets(C);

            //print retained set
            print_retainedset(R);

            printf("\n\n");
        }

        first_round = false;
        // if (data_buffer != NULL){
        //     free(data_buffer);
        // }
    }while(stop_criteria == false);

    // decide what to do with remaining compressed sets and retained set:
    // for now, we will treat them as outliers.
    if(DEBUG){
        //print clusters
        print_clusters(clusters);

        //print compressed sets
        print_compressedsets(C);

        //print retained set
        print_retainedset(R);
    }
    
    if(DEBUG) printf("Freeing data.\n");
    free(clusters);

    
    free(R.points);

    free(C.sets);
    fclose((FILE*)stream_cursor);
    if(DEBUG) printf("Execution complete!.\n");
    return 0;
}