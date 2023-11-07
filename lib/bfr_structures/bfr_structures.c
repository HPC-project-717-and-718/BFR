#include "bfr_structures.h"

RetainedSet init_retained_set(){
    /*
    * Initialize retained set
    *
    * Algorithm:
    *   1. allocate memory for retained set
    *   2. initialize retained set as NULL
    *
    * Parameters:
    *   - void
    *
    * Returns:
    *   - retained set
    */
    RetainedSet R;
    R.points = NULL;
    R.number_of_points = 0;
    return R;
}

CompressedSets init_compressed_sets(){
    /*
    * Initialize compressed sets
    *
    * Algorithm:
    *   1. allocate memory for compressed sets
    *   2. initialize each compressed set as NULL
    *
    * Parameters:
    *   - void
    *
    * Returns:
    *   - array of compressed sets
    */
    CompressedSets C;
    C.sets = NULL;
    C.number_of_sets = 0;
    return C;
}

Cluster * init_cluster(int k){
    /*
    * Initialize clusters
    *
    * Algorithm:
    *   1. allocate memory for clusters
    *   2. initialize each cluster
    *
    * Parameters:
    *   - void
    *
    * Returns:
    *   - array of clusters
    */
    Cluster * clusters = malloc(k * sizeof(Cluster));
    
    if (clusters == NULL){
        printf("Error: could not allocate memory\n");
        exit(1);
    }

    int i = 0;
    for (i = 0; i < k; i++){
        Cluster  c;
        int j = 0;
        for (j = 0; j < M; j++){
            c.centroid.coords[j] = 0.;
            c.sum[j] = 0.;
            c.sum_squares[j] = 0.;
        }
        c.size = 0;
        c.index = i;
        clusters[i] = c;
    }

    return clusters;
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

PriorityQueue* createPriorityQueue(int capacity) {
    /*
    * Create priority queue
    *
    * Algorithm:
    *   1. allocate memory for priority queue
    *   2. allocate memory for priority queue data
    *   3. initialize priority queue size
    *   4. initialize priority queue capacity
    *
    * Parameters:
    *   - capacity: capacity of priority queue
    *
    * Returns:
    *   - pointer to priority queue
    */
    PriorityQueue* pq = (PriorityQueue*)malloc(sizeof(PriorityQueue));
    pq->data = (hierc_element*)malloc(capacity * sizeof(hierc_element));
    pq->size = 0;
    pq->capacity = capacity;
    return pq;
}

void print_clusters(Cluster * clusters){
    printf("Clusters:\n");
    int i = 0;
    for (i = 0; i < K; i++){
        printf("Cluster %d: ", i);
        int j = 0;
        for (j = 0; j < M; j++){
            printf("%lf ", clusters[i].centroid.coords[j]);
        }
        printf("\n");

        printf("Cluster %d size: %d\n", i, clusters[i].size);

        printf("Cluster %d sum: ", i);

        for (j = 0; j < M; j++){
            printf("%lf ", clusters[i].sum[j]);
        }

        printf("\n");

        printf("Cluster %d sum_squares: ", i);

        for (j = 0; j < M; j++){
            printf("%lf ", clusters[i].sum_squares[j]);
        }

        printf("\n");
    }
}

void print_compressedsets(CompressedSets C){
    printf("Compressed Sets:\n");
    int i = 0;
    for(; i < C.number_of_sets; i++){
        printf("Set %d size: %d\n", i, C.sets[i].number_of_points);

        printf("Set %d sum: ", i);
        int j = 0;
        for (j = 0; j < M; j++){
            printf("%d ", C.sets[i].sum[j]);
        }
        printf("\n");

        printf("Set %d sum_squares: ", i);
        j = 0;
        for (j = 0; j < M; j++){
            printf("%d ", C.sets[i].sum_square[j]);
        }
        printf("\n");
    }
}

void print_retainedset(RetainedSet R){
    printf("Retained Set:\n");
    printf("Retained Set size: %d\n", R.number_of_points);
    int i = 0;
    for (; i < R.number_of_points; i++){
        printf("Point %d: ", i);
        int j = 0;
        for (j = 0; j < M; j++){
            printf("%lf ", R.points[i].coords[j]);
        }
        printf("\n");
    }
}


void add_point_to_retained_set(RetainedSet * R, Point p){
    (*R).number_of_points += 1;
    (*R).points = realloc(R->points, R->number_of_points * sizeof(Point));
    if (R->points == NULL){
        printf("Error: could not allocate memory\n");
        exit(1);
    }
    (*R).points[R->number_of_points - 1] = p;
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

void add_miniclusters_to_compressedsets(CompressedSets * C, Cluster * miniclusters, int number_of_miniclusters){
    int i;
    for (i=0; i<number_of_miniclusters; i++){
        if (miniclusters[i].size > 1){
            // add minicluster to compressed sets
            (*C).number_of_sets += 1;
            (*C).sets = realloc(C->sets, C->number_of_sets * sizeof(CompressedSet));
            if (C->sets == NULL){
                printf("Error: could not allocate memory\n");
                exit(1);
            }
            (*C).sets[C->number_of_sets - 1].number_of_points = miniclusters[i].size;
            int j;
            for (j=0; j<M; j++){
                (*C).sets[C->number_of_sets - 1].sum[j] = miniclusters[i].sum[j];
                (*C).sets[C->number_of_sets - 1].sum_square[j] = miniclusters[i].sum_squares[j];
            }
        }
    }
}

CompressedSet merge_compressedsets(CompressedSet C1, CompressedSet C2){
    CompressedSet C;
    C.number_of_points = C1.number_of_points + C2.number_of_points;
    int i;
    for (i=0; i<M; i++){
        C.sum[i] = C1.sum[i] + C2.sum[i];
        C.sum_square[i] = C1.sum_square[i] + C2.sum_square[i];
    }
    return C;
}

void remove_compressedset(CompressedSets * C, int i, int j){
    int k;
    for (k=j; k<C->number_of_sets-1; k++){
        C->sets[k] = C->sets[k+1];
    }
    C->number_of_sets -= 1;
    C->sets = realloc(C->sets, C->number_of_sets * sizeof(CompressedSet));
    if (C->sets == NULL){
        printf("Error: could not allocate memory\n");
        exit(1);
    }

    int l;
    for (l=i; l<C->number_of_sets-1; l++){
        C->sets[l] = C->sets[l+1];
    }
    C->number_of_sets -= 1;
    C->sets = realloc(C->sets, C->number_of_sets * sizeof(CompressedSet));
    if (C->sets == NULL){
        printf("Error: could not allocate memory\n");
        exit(1);
    }
}

void add_compressedset(CompressedSets * C, CompressedSet C1){
    (*C).number_of_sets += 1;
    (*C).sets = realloc(C->sets, C->number_of_sets * sizeof(CompressedSet));
    if (C->sets == NULL){
        printf("Error: could not allocate memory\n");
        exit(1);
    }
    (*C).sets[C->number_of_sets - 1] = C1;
}

bool pop_from_pqueue(PriorityQueue *pq){
    /*
    * Pop from priority queue
    */

    // If the priority queue is empty, return NULL.
    if (pq->size == 0) {
        return false;
    }

    // Save the top of the priority queue.
    hierc_element * hd = malloc(sizeof(hierc_element));
    hd->distance = pq->data[0].distance;
    hd->index_of_cset_1 = pq->data[0].index_of_cset_1;
    hd->index_of_cset_2 = pq->data[0].index_of_cset_2;

    // Move the last element to the top of the priority queue.
    pq->data[0] = pq->data[pq->size - 1];

    // Decrease the size of the priority queue.
    pq->size -= 1;

    // "Sift down" the element that was moved to the top.
    int i = 0;
    while (2 * i + 1 < pq->size) {
        // Find the minimum child.
        int min_child = 2 * i + 1;
        if (2 * i + 2 < pq->size && pq->data[2 * i + 2].distance < pq->data[min_child].distance) {
            min_child = 2 * i + 2;
        }

        // If the parent is greater than the minimum child, swap them.
        if (pq->data[i].distance > pq->data[min_child].distance) {
            hierc_element temp = pq->data[i];
            pq->data[i] = pq->data[min_child];
            pq->data[min_child] = temp;
        }else{
            break;
        }

        // Move to the child.
        i = min_child;
    }

    if (DEBUG) printf("Popped from priority queue: %lf\n", hd->distance);

    return true;
}

bool remove_from_pqueue(PriorityQueue *pq, int i1, int i2){
    PriorityQueue * pq2 = createPriorityQueue(pq->capacity);
    pq2->size = 0;
    int j = 0;
    for (j = 0; j < pq->size; j++){
        if (pq->data[j].index_of_cset_1 != i1 && pq->data[j].index_of_cset_2 != i2 && pq->data[j].index_of_cset_1 != i2 && pq->data[j].index_of_cset_2 != i1){
            add_to_pqueue(pq2, pq->data[j]);
        }
    }
    free(pq->data);
    free(pq);
    
    pq->data = pq2->data;
    pq->size = pq2->size;
}

bool add_to_pqueue(PriorityQueue *pq, hierc_element hd){
    /*
    * Add to priority queue
    *
    * Algorithm:
    *   1. if priority queue is full, return
    *   2. add data to priority queue
    *   3. increase priority queue size
    *   4. sift up priority queue data
    *
    * Parameters:
    *   - pq: priority queue
    *   - hd: data to add to priority queue
    *
    * Returns:
    *   - void
    */
    if (pq->size == pq->capacity){
        return false;
    }
    pq->data[pq->size] = hd;
    pq->size += 1;

    // "Sift up" the new element.
    int i = pq->size - 1;
    while (i != 0 && pq->data[i].distance < pq->data[(i - 1) / 2].distance) {
        // Swap the new element with its parent.
        hierc_element temp = pq->data[i];
        pq->data[i] = pq->data[(i - 1) / 2];
        pq->data[(i - 1) / 2] = temp;

        // Move to the parent.
        i = (i - 1) / 2;
    }
    if (DEBUG) printf("Added to priority queue: %lf\n", hd.distance);
    if (DEBUG){
        //print priority queue in order of priority
        printf("Priority Queue: \n");
        PriorityQueue * pq2 = createPriorityQueue(pq->capacity);
        pq2->size = pq->size;
        int j = 0;
        for (j = 0; j < pq->size; j++){
            pq2->data[j] = pq->data[j];
        }
        while (pq2->size != 0){
            hierc_element * hd2 = top_from_pqueue(*pq2);
            pop_from_pqueue(pq2);
            if ( hd2 == false ){
                printf("error in pop function, the size vlaue do not correspond to real size of vector\n");
                exit(2);
            }
            printf("%lf ", hd2->distance);
        }
        printf("\n");
    }

    return true;
}

bool can_merge(CompressedSet c1, CompressedSet c2){
    //TODO: implement can_merge
    return true;
}

double distance_compressedsets(CompressedSet c1, CompressedSet c2){
    // calculate the two centroids coordinates
    double * coord_c1 = malloc(M * sizeof(double));
    double * coord_c2 = malloc(M * sizeof(double));

    // calculate centroid for c1
    for (int i = 0; i < M; i++) {
        coord_c1[i] = (double)c1.sum[i] / c1.number_of_points;
    }

    // calculate centroid for c2
    for (int i = 0; i < M; i++) {
        coord_c2[i] = (double)c2.sum[i] / c2.number_of_points;
    }

    // calculate the distance between the two centroids
    double distance = 0;
    for (int i = 0; i < M; i++) {
        distance += pow(coord_c1[i] - coord_c2[i], 2);
    }
    distance = sqrt(distance);

    // free the allocated memory
    free(coord_c1);
    free(coord_c2);

    return distance;
}

hierc_element * top_from_pqueue(PriorityQueue pq){
    /*
    * Top from priority queue
    *
    * Algorithm:
    *   1. if priority queue is empty, return NULL
    *   2. return top of priority queue
    *
    * Parameters:
    *   - pq: priority queue
    *
    * Returns:
    *   - top of priority queue
    */
    if (pq.size == 0){
        return NULL;
    }

    hierc_element * hp = malloc(sizeof(hierc_element));
    hp->index_of_cset_1 = pq.data[0].index_of_cset_1;
    hp->index_of_cset_2 = pq.data[0].index_of_cset_2;
    hp->distance = pq.data[0].distance;

    return hp;
}

// unsigned int hash_compressedset(CompressedSet cs) {
//     unsigned int hash = 0;
//     for (int i = 0; i < M; i++) {
//         hash += cs.sum[i];
//         hash += (hash << 10);
//         hash ^= (hash >> 6);
//     }
//     hash += (hash << 3);
//     hash ^= (hash >> 11);
//     hash += (hash << 15);
//     return hash;
// }

bool is_empty_pqueue(PriorityQueue * pq){
    if (pq->size == 0){
        return true;
    }else{
        return false;
    }
}

void merge_compressedsets_and_miniclusters(CompressedSets * C, Cluster * miniclusters, int number_of_miniclusters){
    /*
    * A
    */
    add_miniclusters_to_compressedsets(C, miniclusters, number_of_miniclusters);

    bool stop_merging = false;

    // calculate distances between compressed sets
    PriorityQueue * pq = createPriorityQueue(C->number_of_sets * (C->number_of_sets - 1) / 2);
    int i;
    for (i=0; i<C->number_of_sets; i++){
        int j;
        for (j=i+1; j<C->number_of_sets; j++){
            hierc_element hd;
            hd.distance = distance_compressedsets(C->sets[i], C->sets[j]);
            hd.index_of_cset_1 = i;
            hd.index_of_cset_2 = j;
            add_to_pqueue(pq, hd);
        }
    };

    // naive implementation of hierchical clustering over miniclusters and compressed sets
    while (!stop_merging){
        // top element of priority queue
        hierc_element * hd = top_from_pqueue(*pq);
        if (hd == NULL){
            stop_merging = true;
        }else{
            CompressedSet merged = merge_compressedsets(C->sets[hd->index_of_cset_1], C->sets[hd->index_of_cset_2]);
            pop_from_pqueue(pq);
            if(tightness_evaluation(merged)){
                //TODO: merge and remove from priority queue all element with same index and calculate new distance between merged and all other compressed sets
                // remove_from_pqueue(pq, hd->index_of_cset_1, hd->index_of_cset_2);
            }
        }
        if (is_empty_pqueue(pq)){
            stop_merging = true;
        }
    }

    // TODO: merge compressed sets using hierarchical clustering
    // priority queue pseudocode (naive implementation is also possible, I don't care tbh)
    // priority queue of compressed sets
    // while (priority queue is not empty){
    //     pop two compressed sets
    //     merge them if they fit criteria (to be discussed)
    //     push merged compressed set to priority queue, if it fits criteria
    // }
}