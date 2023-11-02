#include "../lib/kmeans_wrapper/kmeans_wrapper.h"

void take_k_centroids(Cluster *clusters, FILE * input_file){
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
        for (j = 0; j < M; j++){
            
            p.coords[j] = rand() % 10;
        }
        p.cluster = i;
        clusters[i].centroid = p;
    }
}

double mahalanobis_distance(Cluster c, Point p){
    /*
    * Calculates the Mahalanobis radius between a point and a cluster's centroid
    *
    * Algorithm:
    *   1. calculate standard deviation of centroid for each dimension, starting from its variance in said dimension
    *   2. normalize over point coordinate, cluster coordinate and deviation in the dimension
    *   3. square the result, sum over all dimensions and take the root of final result
    *   
    * Parameters:
    *   - c: cluster
    *   - p: point
    *
    * Returns:
    *   - the # of standard deviations of a point from centroid across the dimensions
    */
    double sum = 0, variance = 0, standard_deviation = 0;
    int i = 0;
    for (i = 0; i < M; i++){
        variance = (c.sum_squares[i] / c.size) - pow((c.sum[i] / c.size), 2);
        standard_deviation = sqrt(variance);
        sum += pow(( (p.coords[i] - c.centroid.coords[i]) / standard_deviation ), 2);
    }
    sum = sqrt(sum);

    // printf("distance: %lf\n", sum);
    return sum;
}

double distance(Point p1, Point p2){
    /*
    * Calculate euclidean distance between two points
    *
    * Algorithm:
    *   1. calculate sum of squares of differences between coordinates of two points
    *   2. calculate square root of sum
    *   
    * Parameters:
    *   - p1: first point
    *   - p2: second point
    *
    * Returns:
    *   - distance between two points as double
    */
    double sum = 0;
    int i = 0;
    for (i = 0; i < M; i++){
        sum += (p1.coords[i] - p2.coords[i]) * (p1.coords[i] - p2.coords[i]);
    }
    sum = sqrt(sum);

    // printf("distance: %lf\n", sum);
    return sum;
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

    for (i = 0; i < K; i++){
        current_distance = mahalanobis_distance(clusters[i], p);
        if (current_distance < min_distance) {
            min_distance = current_distance;
            min_cluster = i;
        }
    }

    if (min_distance < T){
        //add point to cluster whose distance from the centroid is minimal
        update_cluster(&clusters[min_cluster], p);
        return true;
    }

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
    double min_distance = distance(clusters[0].centroid, p);
    int index_of_min = 0;
    double second_min_distance = distance(clusters[1].centroid, p);
    int index_of_second_min = 1;
    int i;

    for (i = 1; i < K; i++){
        double dist = distance(clusters[i].centroid, p);
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

    if (distance(pertubated_centroid_1, p) < distance(pertubated_centroid_2, p)){
        //add point to cluster
        update_cluster(&clusters[index_of_min], p);
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

   // should return: cluster with corresponding points, points per cluster
   //cluster_retained_set(R, K);

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
    1. read point from buffer
    2. find if point is close enough to some centroid
    3. if it is, add it to the cluster
    4. add to retained set

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

    printf("streaming data\n");

    Point p;
    int offset = 0;

    while (read_point(data_buffer, &p, size_of_data_buffer, &offset)){
        // apply primary compression criteria
        if (!primary_compression_criteria(clusters, p)){
            // if not, add point to retained set
            (*R).number_of_points += 1;
            (*R).points = realloc(R->points, R->number_of_points * sizeof(Point));
            if (R->points == NULL){
                printf("Error: could not allocate memory\n");
                exit(1);
            }
            (*R).points[R->number_of_points - 1] = p;
        }
    }

    // apply secondary compression criteria over retained set and compressed sets
    secondary_compression_criteria(R, clusters, C);
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

data_streamer data_streamer_Init(char * file_name, char * mode){
    /*
    * Initialize data streamer
    *
    * Algorithm:
    *   1. open file
    *   2. return file pointer
    *
    * Parameters:
    *   - file_name: name of file to open
    *   - mode: mode to open file in
    *
    * Returns:
    *   - file pointer
    */
    FILE * file = fopen(file_name, mode);
    if (file == NULL){
        printf("Error: could not open file\n");
        exit(1);
    }

    //get size of file ---> not needed?
    // fseek(file, 0, SEEK_END);
    // *size_of_file = ftell(file);
    // fseek(file, 0, SEEK_SET);


    return file;
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

    for (i = 0; i < MAX_SIZE_OF_BUFFER; i++){
        Point p;
        for (j = 0; j < M; j++){
            return_val = fscanf(stream_cursor, "%f", &x);

            if (return_val == EOF && i==0 && j==0) {
                // first iteration of the loops, buffer is empty, end algorithm
                return true;
            }
            else if (return_val == EOF) {
                // buffer is not empty - do another iteration, end at the next one.
                return false;
            }

            p.coords[j] = (double) x;
            data_buffer[i] = p;
        }
        counter++;
    }

    *size_of_data_buffer = counter;
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

    data_streamer stream_cursor = data_streamer_Init(argv[1], "r");

    Cluster * clusters = init_cluster(); // allocate the clusters array dynamically

    take_k_centroids(clusters, stream_cursor);

    RetainedSet R = init_retained_set();

    CompressedSets C = init_compressed_sets();

    Point data_buffer[MAX_SIZE_OF_BUFFER];

    //start reading block points from input file
    bool stop_criteria = false;
    do{ 
        long size_of_data_buffer = 0;
        stop_criteria = load_data_buffer(stream_cursor, data_buffer, MAX_SIZE_OF_BUFFER, &size_of_data_buffer);

        if(stop_criteria == false){
            // printf("data buffer: %s\n", data_buffer);
            //stream data
            //TODO: proposition to change the name of this function, in order to make it more clear
            stream_data(&R, clusters, &C, data_buffer, size_of_data_buffer);

            //update centroids
            update_centroids(&clusters, K);

            //print clusters
            print_clusters(clusters);

            //print compressed sets
            print_compressedsets(C);

            //print retained set
            print_retainedset(R);
        }


        // if (data_buffer != NULL){
        //     free(data_buffer);
        // }
    }while(stop_criteria == false);

    // decide what to do with remaining compressed sets and retained set:
    // for now, we will treat them as outliers.
    
    free(clusters);

    
    free(R.points);

    free(C.sets);
    fclose(stream_cursor);
    return 0;
}