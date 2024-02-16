#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include "kmeans.h"


static double pt_distance(const Pointer a, const Pointer b)
{
	Point *pa = (Point*)a;
	Point *pb = (Point*)b;

	double dx = (pa->coords[0] - pb->coords[0]);
	double dy = (pa->coords[1] - pb->coords[1]);

	return dx*dx + dy*dy;
}

static void pt_centroid(const Pointer * objs, const int * clusters, size_t num_objs, int cluster, Pointer centroid)
{
	int i;
	int num_cluster = 0;
	Point sum;
	Point **pts = (Point**)objs;
	Point *center = (Point*)centroid;

	sum.coords[0] = sum.coords[1] = 0.0;

	if (num_objs <= 0) return;

	for (i = 0; i < num_objs; i++)
	{
		/* Only process objects of interest */
		if (clusters[i] != cluster) continue;

		sum.coords[0] += pts[i]->coords[0];
		sum.coords[1] += pts[i]->coords[1];
		num_cluster++;
	}
	if (num_cluster)
	{
		sum.coords[0] /= num_cluster;
		sum.coords[1] /= num_cluster;
		*center = sum;
	}
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
	int i, j;
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
	config.distance_method = pt_distance;
	config.centroid_method = pt_centroid;

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
			double u1 = 1.0 * random() / RAND_MAX;
			double u2 = 1.0 * random() / RAND_MAX;
			double z1 = spread * j + sqrt(-2*log2(u1))*cos(2*M_PI*u2);
			double z2 = spread * j + sqrt(-2*log2(u1))*sin(2*M_PI*u2);
			int n = j*nptsincluster + i;

			/* Populate raw data */
			pts[n].coords[0] = z1;
			pts[n].coords[1] = z2;

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

			if (config.objs[i])
				printf("%g\t%g\t%d\n", pt->coords[0], pt->coords[1], config.clusters[i]);
			else
				printf("N\tN\t%d\n", config.clusters[i]);
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

				if (config.objs[i])
					printf("%g\t%g\t%d\n", pt->coords[0], pt->coords[1], config.clusters[i]);
				else
					printf("N\tN\t%d\n", config.clusters[i]);
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

