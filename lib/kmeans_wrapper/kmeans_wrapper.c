#include "kmeans_wrapper.h"


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
    double sum = 0., variance = 0., standard_deviation = 0.;
    int i = 0;
    for (i = 0; i < M; i++){
        variance = (c.sum_squares[i] / c.size) - pow((c.sum[i] / c.size), 2);
        standard_deviation = sqrt(variance);
        if(standard_deviation < 0.0001) standard_deviation = 1;
        sum += pow(( (p.coords[i] - c.centroid.coords[i]) / standard_deviation ), 2);
        // I was almost going crazy with this part lol
        // if(DEBUG) printf("      Partial2: %lf, Partial1: %lf\n", (c.sum_squares[i] / c.size), pow((c.sum[i] / c.size), 2));
        // if(DEBUG) printf("      SUMSQ: %lf, SUM: %lf, N: %d.\n", c.sum_squares[i], c.sum[i], c.size);
        // if(DEBUG) printf("      P.coords: %lf, Centroid.coords: %lf.\n", p.coords[i], c.centroid.coords[i]);
        // if(DEBUG) printf("      Variance: %lf, Deviation: %lf, sum: %lf.\n", variance, standard_deviation, sum);
        // if(DEBUG) printf("\n");
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
    if(DEBUG) printf("              Updating centroid for cluster %d.\n", cluster);
	int i, j;
	int num_cluster = 0;
	Point sum;
	Point **pts = (Point**)objs;
	Point *center = (Point*)centroid;

	for (i = 0; i < M; i++){
        sum.coords[i] = 0.0;
    }

	if (num_objs <= 0) return;

    if(DEBUG) printf("              Reading all points being clustered.\n");
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
    if(DEBUG) printf("              Read all points, updating cluster coordinates.\n");
	if (num_cluster){
        for (i = 0; i < M; i++){
            sum.coords[i] /= num_cluster;
        }
		*center = sum;
	}
    if(DEBUG) printf("              Cluster coordinates updated.\n");
	return;
}


kmeans_config init_kmeans_config(int k, RetainedSet * R, Point * init, Point * pts){
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

    /* Copy arrays to avoid damaging the RetainedSet*/
	pts = calloc(config.num_objs, sizeof(Point));
	init = calloc(config.k, sizeof(Point));

    if (pts == NULL || init == NULL) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    if(DEBUG) printf("              Copying Retained set.\n");
    /* Add to kmeans' config the RetainedSet's points */
    for (i = 0; i < number_of_points; i++){
        pts[i] = (*R).points[i];
        config.objs[i] = &(pts[i]);
    }
	
    if(DEBUG) printf("              Selecting clusters.\n");
	/* Populate the initial means vector with random start points */
	for (i = 0; i < config.k; i++){
		int r = lround(config.num_objs * (1.0 * rand() / RAND_MAX));
		/* Pointers to the randomly picked point */
        init[i] = (*R).points[r];
		config.centers[i] = &(init[i]);
	}

    return config;
}

Cluster * cluster_retained_set(RetainedSet * R, int k){
    if(DEBUG) printf("          Initializing standard kmeans data.\n");
    Cluster * miniclusters = init_cluster(K);

    Point* init;
    Point* pts;
    kmeans_config config = init_kmeans_config(K, R, init, pts);
    if(DEBUG) printf("          Executing standard kmeans.\n");
	kmeans_result result = kmeans(&config);

    if(DEBUG) printf("          Iteration count: %d\n", config.total_iterations);
    if(DEBUG) printf("          Transferring kmeans cluster data to miniclusters.\n");
    
    // TODO: complete the function by filling miniclusters if a cluster found with kmeans has more than one point
    // TODO: remove from Retained Set the points added to a minicluster, keep the ones not added to a minicluster
    int i;
    for (i = 0; i < config.num_objs; i++){
        Point *pt = (Point *)(config.objs[i]);

        // printf("%g\t%g\t%d\n", pt->coords[0], pt->coords[1], config.clusters[i]);
    }

    if(DEBUG) printf("          Freeing previously allocated data for standard kmeans.\n");
    // free(init); // freeing this data causes a crash. TODO: review error, freeing here probably is not needed as objs[] and clusters[] copy the arrays' data
    // free(pts);  // freeing this data causes a crash. TODO: review error, freeing here probably is not needed as objs[] and clusters[] copy the arrays' data

	free(config.objs);
	free(config.clusters);
	free(config.centers);


    return miniclusters;
}