#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include "kmeans.h"


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
		num_cluster = num_cluster + 1;
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


// the euclidean distance has been modified to fit normal kmeans' requirements
double distance_parallel(const Pointer a, const Pointer b){
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

static void centroid_parallel(const Pointer * objs, const int * clusters, size_t num_objs, int cluster, Pointer centroid){
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
	int num_cluster = 0;
	Point sum;
	Point **pts = (Point**)objs;
	Point *center = (Point*)centroid;

	int i, j;
	// #pragma omp parallel for
	for (i = 0; i < M; i++){
        sum.coords[i] = 0.0;
    }

	if (num_objs <= 0) return;

	// #pragma omp parallel
	// {

	// 	int i, j;
	// 	double sum_coords_private[M];

	// 	#pragma omp parallel for
	// 	for (i = 0; i < M; i++){
	// 		sum_coords_private[i] = 0.0;
	// 	}

	// 	#pragma omp for reduction(+:num_cluster)
	// 	for (i = 0; i < num_objs; i++){
	// 		/* Only process objects of interest */
	// 		if (clusters[i] != cluster) continue;

	// 		for (j = 0; j < M; j++){
	// 			sum_coords_private[j] += pts[i]->coords[j];
	// 		}
	// 		num_cluster = num_cluster + 1;
	// 	}
	// 	#pragma omp critical
	// 	{
	// 		int i;
	// 		for (i = 0; i < M; i++){
	// 			sum.coords[i] += sum_coords_private[i];
	// 		}
	// 	}
	// }

	// #pragma omp parallel for reduction(+:num_cluster,sum_coords[M]) private(j)
	for (i = 0; i < num_objs; i++){
		/* Only process objects of interest */
		if (clusters[i] != cluster) continue;

        // if(DEBUG) printf("              Reading point %d's coordinates.\n", i);
		// #pragma omp parallel for
        for (j = 0; j < M; j++){
            sum.coords[j] += pts[i]->coords[j];
            // if(DEBUG) printf("                  Sum of coords is: %lf\n", sum.coords[i]);
        }
		num_cluster = num_cluster + 1;
	}
    // if(DEBUG) printf("              Read all points, updating cluster coordinates.\n");
	if (num_cluster){
		int i;
		// #pragma omp parallel for
        for (i = 0; i < M; i++){
            sum.coords[i] /= num_cluster;
        }
		*center = sum;
	}
    // if(DEBUG) printf("              Cluster coordinates updated.\n");
	return;
}

int
main(int argc, char **argv)
{
	int rank, size;
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
	kmeans_config config;
	kmeans_result result;
	int i, j, d;
	int spread = 3;
	Point *pts;
	Point *init;
	Point *init_2;
	int print_results = 1;
	unsigned long start;

	int nptsincluster = 100000;
	int k = 10;

	srand(time(NULL));

	/* Constants */
	config.k = k;
	config.num_objs = config.k * nptsincluster;
	config.max_iterations = 200;
	config.distance_method = distance_parallel;
	config.centroid_method = centroid_parallel;

	/* Inputs for K-means */
	config.objs = calloc(config.num_objs, sizeof(Pointer));
	config.centers = calloc(config.k, sizeof(Pointer));
	config.clusters = calloc(config.num_objs, sizeof(int));

	config.parallel = true;
	config.rank = rank;
	config.size = size;

	/* Storage for raw data */
	pts = calloc(config.num_objs, sizeof(Point));
	init = calloc(config.k, sizeof(Point));
	init_2 = calloc(config.k, sizeof(Point));

	/* Create test data! */
	/* Populate with K gaussian clusters of data */
	for (j = 0; j < config.k; j++) {
		for (i = 0; i < nptsincluster; i++)
		{
			int n = j*nptsincluster + i;
			for (d = 0; d < M; d++)
			{
				double u1 = 1.0 * random() / RAND_MAX;
				double u2 = 1.0 * random() / RAND_MAX;
				double z1 = spread * j + sqrt(-2*log2(u1))*cos(2*M_PI*u2);

				/* Populate raw data */
				pts[n].coords[d] = z1;

			}
			/* Pointer to raw data */
			config.objs[n] = &(pts[n]);			
		}
	}

	/* Populate the initial means vector with random start points */
	for (i = 0; i < config.k; i++)
	{
		int r = lround(config.num_objs * (1.0 * rand() / RAND_MAX));
		/* Populate raw data */
		init[i] = pts[r];
		init_2[i] = pts[r];
		/* Pointers to raw data */
		config.centers[i] = &(init[i]);

		Point* pointArray = (Point*)config.centers[i];
		double* coords_ptr = pointArray->coords;

		if (print_results && rank == 0)
			printf("center[%d]\t%g\t%g\n", i, coords_ptr[0], coords_ptr[1]);
	}

	/* run k-means! */
	start = time(NULL);
	result = kmeans(&config);

	/* print results */
	if (print_results && rank == MASTER)
	{

		printf("\n");
		printf("Iteration count: %d\n", config.total_iterations);
		printf("     Time taken: %ld seconds\n", (time(NULL) - start));
		printf(" Iterations/sec: %.3g\n", (1.0*config.total_iterations)/(time(NULL) - start));
		printf("\n");
		for (i = 0; i < config.num_objs; i++)
		{
			Point *pt = (Point*)(config.objs[i]);


			if (config.objs[i]){
				for (d = 0; d < M; d++)
				{
					printf("%g\t", pt->coords[d]);
				}
				printf("%d\n", config.clusters[i]);
			}
			else{
				printf("N\tN\t%d\n", config.clusters[i]);
			}
			
		}
	}

	if(rank == MASTER){
		printf("AGAIN!!!");
		/* Populate the initial means vector with random start points */
		for (i = 0; i < config.k; i++)
		{
			/* Pointers to raw data */
			config.centers[i] = &(init_2[i]);

			Point* pointArray = (Point*)config.centers[i];
			double* coords_ptr = pointArray->coords;

			if (print_results && rank == 0)
				printf("center[%d]\t%g\t%g\n", i, coords_ptr[0], coords_ptr[1]);
		}

		config.parallel = false;
		config.distance_method = distance;
		config.centroid_method = centroid;
		
		/* run k-means! */
		start = time(NULL);
		result = kmeans(&config);

		/* print results */
		if (print_results && rank == MASTER)
		{

			printf("\n");
			printf("Iteration count: %d\n", config.total_iterations);
			printf("     Time taken: %ld seconds\n", (time(NULL) - start));
			printf(" Iterations/sec: %.3g\n", (1.0*config.total_iterations)/(time(NULL) - start));
			printf("\n");
			for (i = 0; i < config.num_objs; i++)
			{
				Point *pt = (Point*)(config.objs[i]);


				if (config.objs[i]){
					for (d = 0; d < M; d++)
					{
						printf("%g\t", pt->coords[d]);
					}
					printf("%d\n", config.clusters[i]);
				}
				else{
					printf("N\tN\t%d\n", config.clusters[i]);
				}
				
			}
		}
	}

	free(config.objs);
	free(config.clusters);
	free(config.centers);

	free(init);
	free(init_2);
	free(pts);
	
    // Finalize the MPI environment
    MPI_Finalize();
    return 0;
}

