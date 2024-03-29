#include "bfr_helper_functions.h"


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
    double sum = 0., sum_partial = 0., variance = 0.;
    int i = 0;
    for (i = 0; i < M; i++){
        variance = (c.sum_squares[i] / c.size) - pow((c.sum[i] / c.size), 2);
        if(variance < 0.0001) variance = 0.0001;
        sum_partial = pow((p.coords[i] - c.centroid.coords[i]) , 2);
        sum = sum + (sum_partial / variance);
        // sum += pow(( (p.coords[i] - c.centroid.coords[i]) / standard_deviation ), 2);
        // I was almost going crazy with this part lol
        // if(DEBUG) printf("\n");
        // if(DEBUG) printf("      Partial2: %lf, Partial1: %lf\n", (c.sum_squares[i] / c.size), pow((c.sum[i] / c.size), 2));
        // if(DEBUG) printf("      SUMSQ: %lf, SUM: %lf, N: %d.\n", c.sum_squares[i], c.sum[i], c.size);
        // if(DEBUG) printf("      P.coords: %lf.\n", p.coords[i]);
        // if(DEBUG) printf("      Centroid.coords: %lf.\n", c.centroid.coords[i]);
        // if(DEBUG) printf("      Variance: %lf, Deviation: %lf, sum: %lf.\n", variance, standard_deviation, sum);
    }
    sum = sqrt(sum);

    // printf("distance: %lf\n", sum);
    return sum;
}

// the euclidean distance has been modified to fit normal kmeans' requirements
double distance(const Pointer a, const Pointer b){
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
    Point* p1 = (Point*)a;
    Point* p2 = (Point*)b;
    //TODO: insert control if the cast is successful or not

    double sum = 0;
    int i = 0;
    for (i = 0; i < M; i++){
        sum += (p1->coords[i] - p2->coords[i]) * (p1->coords[i] - p2->coords[i]);
    }
    sum = sqrt(sum);

    return sum;
}

static void centroid(const Pointer * objs, const int * clusters, size_t num_objs, int cluster, Pointer centroid){
    /*
    * Given a collection of objects categorized by cluster, returns their centroid
    *
    * Algorithm:
    *   1. calculate sum of coordinates in all dimentions for all points of the cluster
    *   2. update cluster coordinates
    *   
    * Parameters:
    *   - objs: all points being clustered
    *   - clusters: array of cluster ids
    *   - num_objs: # of objs
    *   - cluster: cluster that needs update
    *   - centroid: centroid of cluster that needs update
    *
    * Returns:
    *   - void
    * 
    * Note: this might be inefficient, as in BFR we could save a cluster by its statistics directly
    *       but for now we'll stick to normal kmeans' requirements
    */
    // if(DEBUG) printf("              Updating centroid for cluster %d.\n", cluster);
	int i, j;
	int num_cluster = 0;
	Point sum;
	Point **pts = (Point**)objs;
	Point *center = (Point*)centroid;

	for (i = 0; i < M; i++){
        sum.coords[i] = 0.0;
    }

	if (num_objs <= 0) return;

    // if(DEBUG) printf("              Reading all points being clustered.\n");
	for (i = 0; i < num_objs; i++){
		/* Only process objects of interest */
		if (clusters[i] != cluster) continue;

        // if(DEBUG) printf("              Reading point %d's coordinates.\n", i);
        for (j = 0; j < M; j++){
            sum.coords[j] += pts[i]->coords[j];
            // if(DEBUG) printf("                  Sum of coords is: %lf\n", sum.coords[i]);
        }
		num_cluster++;
	}
    // if(DEBUG) printf("              Read all points, updating cluster coordinates.\n");
	if (num_cluster){
        for (i = 0; i < M; i++){
            sum.coords[i] /= num_cluster;
        }
		*center = sum;
	}
    // if(DEBUG) printf("              Cluster coordinates updated.\n");
	return;
}


kmeans_config init_kmeans_config(int k, RetainedSet * R, bool parallel, int rank, int size){
    kmeans_config config;

    int i, number_of_points = (*R).number_of_points;

    if(DEBUG) printf("              Declaring constants.\n");
	/* Constants */
	config.k = k;
	config.num_objs = number_of_points;
	config.max_iterations = NORMAL_KMEANS_MAX_ITERATIONS;
	config.distance_method = distance;
	config.centroid_method = centroid;

    if(DEBUG) printf("              Allocating dynamic data.\n");
	/* Inputs for K-means */
	config.objs = calloc(config.num_objs, sizeof(Pointer));
	config.centers = calloc(config.k, sizeof(Pointer));
	config.clusters = calloc(config.num_objs, sizeof(int));

    config.parallel = parallel;
    config.rank = rank;
    config.size = size;

    if(DEBUG) printf("              Copying Retained set.\n");
    /* Add to kmeans' config the RetainedSet's points */
    for (i = 0; i < number_of_points; i++){

        config.objs[i] = &((*R).points[i]);
    }
	
    if(DEBUG) printf("              Selecting clusters.\n");
	/* Populate the initial means vector with random start points */
	for (i = 0; i < config.k; i++){
		// int r = lround((config.num_objs-1) * (1.0 * rand() / RAND_MAX));
        // if(DEBUG) printf("              r:%d\n", r);
        // if(DEBUG) printf("              R.points[r]:%lf %lf.\n", (*R).points[r].coords[0], (*R).points[r].coords[1]);
		/* Pointers to the randomly picked point */
        
		// config.centers[i] = &((*R).points[r]);
		config.centers[i] = &((*R).points[i]);  // just set centers to i-th point of R, we need the algorithm to be deterministic
	}

    return config;
}

bool tightness_evaluation_cluster(Cluster c, int * tightness_flag, int index){
    /*
    */
    //if tighteness_flag[index] 0 is not computed, compute it
    //if tighteness_flag[index] 1 is computed, and is false
    //if tighteness_flag[index] 2 is computed, and is true
    // if(DEBUG) printf("              Evaluating tightness for minicluster %d.\n", index);
    if (tightness_flag[index] == 0){
        if (c.size < 2){
            tightness_flag[index] = 1;
            return false;
        }
        int x_sub[M];
        int i = 0;
        for (i = 0; i < M; i++){
            x_sub[i] = c.sum[i];
            x_sub[i] = x_sub[i] / c.size;
        }
        int j = 0;
        double max_value = 0;
        for (j = 0; j < M; j++){
            double value = c.sum_squares[j] - (c.size * pow(x_sub[j], 2));
            value = value / (c.size);
            value = sqrt(value);
            if (value > max_value){
                max_value = value;
            }
        }
        // if(DEBUG) printf("              Max value for tightness constraint: %lf.\n", max_value);
        if (max_value < BETA){
            tightness_flag[index] = 2;
            return true;
        }else{
            tightness_flag[index] = 1;
            return false;
        }
    }
    else if (tightness_flag[index] == 1){
        return false;
    }
    else if (tightness_flag[index] == 2){
        return true;
    }
    else {
        printf("Error: tightness_flag[%d] has an invalid value: %d\n", index, tightness_flag[index]);
        exit(1);
    } 

}

int update_miniclusters(Cluster ** miniclusters, int * tightness_flag, int k){
    /*
    */
    int i = 0;
    int j = 0;
    int new_k = 0;
    for (i = 0; i < k; i++){
        if (tightness_flag[i] == 2){
            new_k++;
        }
    }
    Cluster * new_miniclusters = init_cluster(new_k);
    j = 0;
    for (i = 0; i < k; i++){
        if (tightness_flag[i] == 2){
            new_miniclusters[j] = (*miniclusters)[i];
            j++;
        }
    }
    free(*miniclusters);
    *miniclusters = new_miniclusters;

    return new_k;
}

Cluster * cluster_retained_set(RetainedSet * R, int * k){
    if(DEBUG) printf("          Initializing standard kmeans data.\n");
    Cluster * miniclusters = init_cluster((*k));

    // TODO: discuss a correct limit for not running standard kmeans
    // as of now, do not run kmeans if the number of points is < k
    // we may want to run kmeans when we have more than k*constant number of points
    if((*R).number_of_points < (*k)){
        if(DEBUG) printf("          Retained set has less points (%d) than clusters(%d). Returning empty miniclusters.\n", (*R).number_of_points, (*k));
        return miniclusters;
    }

    kmeans_config config = init_kmeans_config((*k), R, false, 0, 0);
    if(DEBUG) printf("          Executing standard kmeans.\n");
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
    
    return miniclusters;
}

bool load_data_buffer(data_streamer stream_cursor, Point * data_buffer, int max_size_of_buffer, long * size_of_data_buffer){
    if (stream_cursor == NULL){
        printf("Error: invalid stream cursor\n");
        exit(1);
    }
    
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
                *size_of_data_buffer = counter;
                if(DEBUG) printf("  Data buffer size is %d.\n", counter);
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

    return false;
}