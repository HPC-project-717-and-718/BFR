/*-------------------------------------------------------------------------
*
* kmeans.c
*    Generic k-means implementation
*
* Copyright (c) 2016, Paul Ramsey <pramsey@cleverelephant.ca>
*
*------------------------------------------------------------------------*/


#include "kmeans.h"

static void
update_r(kmeans_config *config)
{
	int i;

	for (i = 0; i < config->num_objs; i++)
	{
		double distance, curr_distance;
		int cluster, curr_cluster;
		Pointer obj;

		assert(config->objs != NULL);
		assert(config->num_objs > 0);
		assert(config->centers);
		assert(config->clusters);

		obj = config->objs[i];

		/*
		* Don't try to cluster NULL objects, just add them
		* to the "unclusterable cluster"
		*/
		if (!obj)
		{
			config->clusters[i] = KMEANS_NULL_CLUSTER;
			continue;
		}

		/* Initialize with distance to first cluster */
		curr_distance = (config->distance_method)(obj, config->centers[0]);
		curr_cluster = 0;

		/* Check all other cluster centers and find the nearest */
		for (cluster = 1; cluster < config->k; cluster++)
		{
			distance = (config->distance_method)(obj, config->centers[cluster]);
			if (distance < curr_distance)
			{
				curr_distance = distance;
				curr_cluster = cluster;
			}
		}

		/* Store the nearest cluster this object is in */
		config->clusters[i] = curr_cluster;
	}
}

static void
update_means(kmeans_config *config)
{
	int i;

	for (i = 0; i < config->k; i++)
	{
		/* Update the centroid for this cluster */
		(config->centroid_method)(config->objs, config->clusters, config->num_objs, i, config->centers[i]);
	}
}

#ifdef KMEANS_THREADED

// static void * update_r_threaded_main(void *args)
// {
// 	kmeans_config *config = (kmeans_config*)args;
// 	update_r(config);
// 	pthread_exit(args);
// }

// static void update_r_threaded(kmeans_config *config)
// {
// 	/* Computational complexity is function of objs/clusters */
// 	/* We only spin up threading infra if we need more than one core */
// 	/* running. We keep the threshold high so the overhead of */
// 	/* thread management is small compared to thread compute time */
// 	int num_threads = config->num_objs * config->k / KMEANS_THR_THRESHOLD;

// 	/* Can't run more threads than the maximum */
// 	num_threads = (num_threads > KMEANS_THR_MAX ? KMEANS_THR_MAX : num_threads);

// 	/* If the problem size is small, don't bother w/ threading */
// 	if (num_threads < 1)
// 	{
// 		update_r(config);
// 	}
// 	else
// 	{
// 		pthread_t thread[KMEANS_THR_MAX];
// 		pthread_attr_t thread_attr;
// 		kmeans_config thread_config[KMEANS_THR_MAX];
// 		int obs_per_thread = config->num_objs / num_threads;
// 		int i, rc;

// 		for (i = 0; i < num_threads; i++)
// 		{
// 			/*
// 			* Each thread gets a copy of the config, but with the list pointers
// 			* offest to the start of the batch the thread is responsible for, and the
// 			* object count number adjusted similarly.
// 			*/
// 			memcpy(&(thread_config[i]), config, sizeof(kmeans_config));
// 			thread_config[i].objs += i*obs_per_thread;
// 			thread_config[i].clusters += i*obs_per_thread;
// 			thread_config[i].num_objs = obs_per_thread;
// 			if (i == num_threads-1)
// 			{
// 				thread_config[i].num_objs += config->num_objs - num_threads*obs_per_thread;
// 			}

// 			/* Initialize and set thread detached attribute */
// 			pthread_attr_init(&thread_attr);
// 			pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

// 			/* Now we just run the thread, on its subset of the data */
// 			rc = pthread_create(&thread[i], &thread_attr, update_r_threaded_main, (void *) &thread_config[i]);
// 			if (rc)
// 			{
				// printf("ERROR: return code from pthread_create() is %d\n", rc);
// 				exit(-1);
// 			}
// 		}

// 		/* Free attribute and wait for the other threads */
// 		pthread_attr_destroy(&thread_attr);

// 		/* Wait for all calculations to complete */
// 		for (i = 0; i < num_threads; i++)
// 		{
// 		    void *status;
// 			rc = pthread_join(thread[i], &status);
// 			if (rc)
// 			{
				// printf("ERROR: return code from pthread_join() is %d\n", rc);
// 				exit(-1);
// 			}
// 		}
// 	}
// }

static void * update_r_parallel_main(void *args)
{
	kmeans_config *config = (kmeans_config*)args;
	update_r(config);
}

static void update_r_parallel(kmeans_config *config)
{
	int obs_per_node = config->num_objs / config->size;

	/* For each node, create a config copy with the objects, clusters and num_objs offset correctly*/
	kmeans_config node_config;

	// if(config->rank == MASTER){
	// 	// printf("0: cluster assignment BEFORE\n");
	// 	int i;
	// 	for (i = 0; i < config->num_objs; i++){
	// 		// printf("%d ", config->clusters[i]);
	// 	}
	// 	// printf("\n");
	// }

	// printf("%d: Copying memory.\n", config->rank);
	memcpy(&(node_config), config, sizeof(kmeans_config));
	node_config.objs += config->rank*obs_per_node;
	node_config.clusters += config->rank*obs_per_node;
	node_config.num_objs = obs_per_node;
	if (config->rank == config->size-1)
	{
		node_config.num_objs += config->num_objs - config->size*obs_per_node;
	}

	// printf("%d: Updating r main.\n", config->rank);
	/* Run the node, on its subset of the data */
	update_r_parallel_main((void *) &node_config);

	// Determine displacement and count for each process
    int* displs = (int*)malloc(sizeof(int) * config->size);
    int* rcounts = (int*)malloc(sizeof(int) * config->size);
    int total_elements = 0, i;
	for (i = 0; i < config->size; ++i) {
        displs[i] = total_elements;
        rcounts[i] = obs_per_node;
		if (i == config->size-1){
			rcounts[i] += config->num_objs - config->size*obs_per_node;
		}
        total_elements += node_config.num_objs;
    }

	int *clusters_copy = malloc(node_config.num_objs*sizeof(int));

	// printf("%d: Printing local results.\n",config->rank);
	for (i = 0; i < node_config.num_objs; i++) {
		// printf("%d: %d\n", config->rank, node_config.clusters[i]);
		clusters_copy[i] = node_config.clusters[i];
	}
	// printf("\n\n");

	// printf("%d: Copying local results.\n",config->rank);
	// memcpy(clusters_copy, node_config.clusters, rcounts[config->rank]*sizeof(int));
	
	// printf("%d: Gathering results.\n", config->rank);
	// MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, config->clusters, rcounts, displs, MPI_INT, MPI_COMM_WORLD);
	MPI_Allgatherv(clusters_copy, rcounts[config->rank], MPI_INT, config->clusters, rcounts, displs, MPI_INT, MPI_COMM_WORLD);	

	// if(config->rank == MASTER){
	// // printf("%d: cluster assignment AFTER\n", config->rank);
	// for (i = 0; i < config->num_objs; i++){
	// 	// printf("%d ", config->clusters[i]);
	// }
	// printf("\n");
	// }

	free(clusters_copy);
}

// int update_means_k;
// pthread_mutex_t update_means_k_mutex;

// static void *
// update_means_threaded_main(void *arg)
// {
// 	kmeans_config *config = (kmeans_config*)arg;
// 	int i = 0;

// 	do
// 	{
// 		pthread_mutex_lock (&update_means_k_mutex);
// 		i = update_means_k;
// 		update_means_k++;
// 		pthread_mutex_unlock (&update_means_k_mutex);

// 		if (i < config->k)
// 			(config->centroid_method)(config->objs, config->clusters, config->num_objs, i, config->centers[i]);
// 	}
// 	while (i < config->k);

// 	pthread_exit(arg);
// }

// static void
// update_means_threaded(kmeans_config *config)
// {
// 	/* We only spin up threading infra if we need more than one core */
// 	/* running. We keep the threshold high so the overhead of */
// 	/* thread management is small compared to thread compute time */
// 	int num_threads = config->num_objs / KMEANS_THR_THRESHOLD;

// 	/* Can't run more threads than the maximum */
// 	num_threads = (num_threads > KMEANS_THR_MAX ? KMEANS_THR_MAX : num_threads);

// 	/* If the problem size is small, don't bother w/ threading */
// 	if (num_threads < 1)
// 	{
// 		update_means(config);
// 	}
// 	else
// 	{
// 		/* Mutex protected counter to drive threads */
// 		pthread_t thread[KMEANS_THR_MAX];
// 		pthread_attr_t thread_attr;
// 		int i, rc;

// 		pthread_mutex_init(&update_means_k_mutex, NULL);
// 		update_means_k = 0;

// 		pthread_attr_init(&thread_attr);
// 		pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

// 		/* Create threads to perform computation  */
// 		for (i = 0; i < num_threads; i++)
// 		{

// 			/* Now we just run the thread, on its subset of the data */
// 			rc = pthread_create(&thread[i], &thread_attr, update_means_threaded_main, (void *) config);
// 			if (rc)
// 			{
				// printf("ERROR: return code from pthread_create() is %d\n", rc);
// 				exit(-1);
// 			}
// 		}

// 		pthread_attr_destroy(&thread_attr);

// 		/* Watch until completion  */
// 		for (i = 0; i < num_threads; i++)
// 		{
// 		    void *status;
// 			rc = pthread_join(thread[i], &status);
// 			if (rc)
// 			{
				// printf("ERROR: return code from pthread_join() is %d\n", rc);
// 				exit(-1);
// 			}
// 		}

// 		pthread_mutex_destroy(&update_means_k_mutex);
// 	}
// }

static void
update_means_parallel(kmeans_config *config)
{
	/* What we want to do here is assign each node a number of clusters. Each will update its center and the results will be shared to all nodes. */
	if(config->k <= config->size){
		// printf("%d: Means are less than nodes.\n", config->rank);

		if(config->rank < config->k){
			// printf("%d: Receiving one mean to edit.\n", config->rank);
			/* Assign one cluster to each node until no cluster can be assigned */
			int clusters_per_node = 1;
			int clusters_per_node_all = clusters_per_node;

			/* For each node, recompute mean and send new centroid coordinates to MASTER */
			int offset;
			offset = config->rank;
			// printf("%d: Calling centroid method for mean of offset %d. Coords were %lf, %lf.\n", config->rank, offset, ((Point*)config->centers[offset])->coords[0], ((Point*)config->centers[offset])->coords[1]);
			(config->centroid_method)(config->objs, config->clusters, config->num_objs, offset, config->centers[offset]);
			
			Point* pointArray = (Point*)config->centers[offset];
	
			double* coords_ptr = pointArray->coords;

			// printf("%d: Point coords are %lf, %lf.\n", config->rank, coords_ptr[0], coords_ptr[1]);
			if(config->rank != MASTER) {
				// printf("%d: Sending coords to MASTER.\n", config->rank);
				MPI_Send(coords_ptr, M, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);
			}
		}

		/* MASTER has to receive each centroid coordinates and update them accordingly */
		if(config->rank == MASTER){
			// printf("%d: Receiving coords (MASTER).\n", config->rank);
			double* temp_coords_ptr = (double*)malloc(M*sizeof(double));
			/* Receive new centers */
			int i, offset;
			for(i = 1; i < config->k; i++){
				offset = config->rank;
				MPI_Recv(temp_coords_ptr, M, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
				Point* pointArray = (Point*)config->centers[offset];
				// Copy received values to the coords array
				int k;
				for (k = 0; k < M; k++) {
					pointArray->coords[k] = temp_coords_ptr[k];
    			}
			}
			// printf("%d: Freeing temp coords ptr.\n", config->rank);
			free(temp_coords_ptr);
		}

	}
	else{
		// printf("%d: There are more means than nodes.\n", config->rank);
		/* Assign a (possibly) equal number of clusters to each node*/
		int clusters_per_node = config->k / config->size;
		int clusters_per_node_all = clusters_per_node;
		if (config->rank == config->size-1)
		{
			/* Last node gets the remainder */
			clusters_per_node += config->k - config->size*clusters_per_node;
		}

		/* For each node, recompute mean and send new centroid coordinates to MASTER */
		int i, offset;
		for(i = 0; i < clusters_per_node; i++){
			offset = i + config->rank * clusters_per_node_all;
			(config->centroid_method)(config->objs, config->clusters, config->num_objs, offset, config->centers[offset]);
            
			Point* pointArray = (Point*)config->centers[offset];
	
			double* coords_ptr = pointArray->coords;
			if(config->rank != MASTER) MPI_Send(coords_ptr, M, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD);
		}

		/* MASTER has to receive each centroid coordinates and update them accordingly */
		if(config->rank == MASTER){
			double* temp_coords_ptr = (double*)malloc(M*sizeof(double));
			/* Receive new centers */
			int j;
			for(i = 1; i < config->size; i++){
				clusters_per_node = config->k / config->size;
				clusters_per_node_all = clusters_per_node;
				if (config->rank == config->size-1){
					clusters_per_node += config->k - config->size*clusters_per_node;
				}

				for(j = 0; j < clusters_per_node; j++){
					offset = j + config->rank * clusters_per_node_all;
					MPI_Recv(temp_coords_ptr, M, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					
					Point* pointArray = (Point*)config->centers[offset];

					// Copy received values to the coords array
					int k;
					for (k = 0; k < M; k++) {
						pointArray->coords[k] = temp_coords_ptr[k];
    				}
				}
			}
			free(temp_coords_ptr);
		}
	}

	// printf("%d: Broadcasting coordinates.\n", config->rank);
	int i, offset;
	for (i = 0; i < config->k; i++){
		Point* pointArray = (Point*)config->centers[i];
		double* coords_ptr = pointArray->coords;
		MPI_Bcast(coords_ptr, M, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
		if(config->rank != MASTER){
			// printf("%d: Updating coordinates.\n", config->rank);
			/* Update centroid coordinates */
			int k;
			for (k = 0; k < M; k++) {
				pointArray->coords[k] = coords_ptr[k];
    		}
		}
	}
}

#endif /* KMEANS_THREADED */

kmeans_result
kmeans(kmeans_config *config)
{
	int iterations = 0;
	int *clusters_last;
	size_t clusters_sz = sizeof(int)*config->num_objs;

	assert(config);
	assert(config->objs);
	assert(config->num_objs);
	assert(config->distance_method);
	assert(config->centroid_method);
	assert(config->centers);
	assert(config->k);
	assert(config->clusters);
	assert(config->k <= config->num_objs);
	// assert(config->parallel);
	// assert(config->rank);
	// assert(config->size);

	/* Zero out cluster numbers, just in case user forgets */
	memset(config->clusters, 0, clusters_sz);
	// printf("\n%d: Starting kmeans.\n", config->rank);

	/* Set default max iterations if necessary */
	if (!config->max_iterations)
		config->max_iterations = KMEANS_MAX_ITERATIONS;

	/*
	 * Previous cluster state array. At this time, r doesn't mean anything
	 * but it's ok
	 */
	clusters_last = kmeans_malloc(clusters_sz);

	while (1)
	{
		if(config->rank == MASTER){
			printf("0: cluster assignment BEFORE\n");
			int i;
			for (i = 0; i < config->num_objs; i++){
				printf("%d ", config->clusters[i]);
			}
			printf("\n");
		}
		// printf("\n%d: Iterating.\n", config->rank);
		/* Store the previous state of the clustering */
		memcpy(clusters_last, config->clusters, clusters_sz);

#ifdef KMEANS_THREADED
		if(config->parallel){
			// printf("\n%d: In parallel, updating r.\n", config->rank);
			/* At this point, all nodes have the same config. Have master coordinate the clustering, then broadcast the results. */
			update_r_parallel(config);
			// printf("\n%d: In parallel, updating means.\n", config->rank);
			update_means_parallel(config);
		}
		else{
#endif
			update_r(config);
			update_means(config);
#ifdef KMEANS_THREADED
		}
#endif
		if(config->rank == MASTER){
			printf("0: cluster assignment AFTER\n");
			int i;
			for (i = 0; i < config->num_objs; i++){
				printf("%d ", config->clusters[i]);
			}
			printf("\n");
		}
		/*
		 * if all the cluster numbers are unchanged since last time,
		 * we are at a stable solution, so we can stop here
		 */
		if (memcmp(clusters_last, config->clusters, clusters_sz) == 0)
		{
			kmeans_free(clusters_last);
			config->total_iterations = iterations;
			return KMEANS_OK;
		}

		if (iterations++ > config->max_iterations)
		{
			kmeans_free(clusters_last);
			config->total_iterations = iterations;
			return KMEANS_EXCEEDED_MAX_ITERATIONS;
		}
	}

	kmeans_free(clusters_last);
	config->total_iterations = iterations;
	return KMEANS_ERROR;
}


