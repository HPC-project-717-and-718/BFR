#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#define K 2 // number of clusters
#define M 2 // number of dimensions
#define S 100 // sample size
#define T 3.5 // threshold
#define LIMIT_S 10 //limit for rand generation of coordinates for initial centroids

typedef struct {
    double coords[M];
    int cluster;
} Point;

typedef struct {
    Point centroid;
    int size;
    int sum[M];
    int sum_squares[M];
    int index;
} Cluster;

typedef struct {
    Point * points;
    int number_of_points;
} RetainedSet;

typedef struct {
    int number_of_points;
    int sum[M];
    int sum_square[M];
} CompressedSet;

typedef struct {
    CompressedSet * sets;
    int number_of_sets;
} CompressedSets;

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

Cluster * init_cluster(){
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
    Cluster * clusters = malloc(K * sizeof(Cluster));
    
    if (clusters == NULL){
        printf("Error: could not allocate memory\n");
        exit(1);
    }

    int i = 0;
    for (i = 0; i < K; i++){
        Cluster  c;
        int j = 0;
        for (j = 0; j < M; j++){
            c.centroid.coords[j] = 0;
            c.sum[j] = 0;
            c.sum_squares[j] = 0;
        }
        c.size = 0;
        c.index = i;
        clusters[i] = c;
    }

    return clusters;
}

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

    printf("distance: %lf\n", sum);
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

bool read_point(FILE * input_file, Point * p){
    /*
    * Read a point from input file
    *
    * Algorithm:
    *   1. read M coordinates from input file
    *   2. if EOF is reached return false
    *   3. else return true
    *
    * Parameters:
    *   - input_file: file to read point from
    *   - p: point to read
    *
    * Returns:
    *   - true if point is read successfully
    *   - false if EOF is reached
    */
    int i = 0;
    for (i = 0; i < M; i++){
        if (fscanf(input_file, "%lf", &p->coords[i]) == EOF){
            return false;
        }
    }
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

void stream_data(RetainedSet * R, Cluster * clusters, CompressedSets * C, FILE * input_file){

    /*
    1. read point from input file
    2. find if point is close enough to some centroid
    3. if it is, add it to the cluster
    4. if it is not, decide if add point to compressed set or retained set

    * Parameters:
    *   - R: retained set
    *   - clusters: array of clusters
    *   - C: compressed sets
    *   - input_file: file to read points from
    *
    * Returns:
    *   - void
    */

    printf("streaming data\n");
    Point p;
    while (read_point(input_file, &p)){
        // find if point is close enough to some centroid
        // if it is, add it to the cluster
        bool flag = false;
        int i = 0;
        for (i = 0; i < K && !flag; i++){
            if (distance(p, clusters[i].centroid) < T){
                printf("point is close to centroid\n");
                flag = true;

                update_cluster(&clusters[i], p);
            }
        }

        if (!flag){
            // decide if add point to compressed set or retained set
            // iter over compressed set decide if add to one of the sets if distance between the center of the set and the point is less than T
            // the center of the set is the coordinates in each dimension divided by number of points in the set
            // if point is not added to any set, add it to retained set
            printf("point is not close to any centroid\n");

            if (C->number_of_sets > 0){
                int size_of_C = C->number_of_sets;
                int j = 0;
                for (j = 0; j < size_of_C; j++){
                    //calculate center of set
                    Point center;
                    int k = 0;
                    for (k = 0; k < M; k++){
                        center.coords[k] = C->sets[j].sum[k] / C->sets[j].number_of_points;
                    }
                    if (distance(p, center) < T){
                        //add point to set
                        C->sets[j].number_of_points += 1;
                        int k = 0;
                        for (k = 0; k < M; k++){
                            C->sets[j].sum[k] += p.coords[k];
                            C->sets[j].sum_square[k] += pow(p.coords[k], 2);
                        }
                        break;
                    }
                }
            }else{
                //check if point is close to any point in retained set
                // in this case merge the two points and create a new compressed set
                bool is_close = false;
                int size_of_R = R->number_of_points;
                int j = 0;
                for (j = 0; j < size_of_R; j++){
                    if (distance(p, R->points[j]) < T){
                        //merge points
                        merge_points(&p, &R->points[j], C, R);
                        is_close = true;
                        break;
                    }
                }

                if(!is_close){
                    //add point to retained set
                    (*R).number_of_points += 1;
                    (*R).points = realloc(R->points, R->number_of_points * sizeof(Point));
                    if (R->points == NULL){
                        printf("Error: could not allocate memory\n");
                        exit(1);
                    }
                    (*R).points[R->number_of_points - 1] = p;
                }
                
            }

        }
    }
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
            printf("%d ", clusters[i].sum[j]);
        }

        printf("\n");

        printf("Cluster %d sum_squares: ", i);

        for (j = 0; j < M; j++){
            printf("%d ", clusters[i].sum_squares[j]);
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

int main(int argc, char ** argv){
    if (argc != 2){
        printf("Usage: %s <input_file>\n", argv[0]);
        return 1;
    }

    FILE * input_file = fopen(argv[1], "r");

    if (input_file == NULL){
        printf("Error opening file %s\n", argv[1]);
        return 1;
    }

    Cluster * clusters = init_cluster(); // allocate the clusters array dynamically

    take_k_centroids(clusters, input_file);

    printf("Centroids:\n");
    int i = 0;
    for (i = 0; i < K; i++){
        printf("Centroid %d: ", i);
        int j = 0;
        for (j = 0; j < M; j++){
            printf("%lf ", clusters[i].centroid.coords[j]);
        }
        printf("\n");
    }

    RetainedSet R = init_retained_set();

    CompressedSets C = init_compressed_sets();
    //start reading points from input file
    stream_data(&R, clusters, &C, input_file);

    print_clusters(clusters);
    print_compressedsets(C);
    print_retainedset(R);
     //update centroids
    update_centroids(&clusters, K);

    print_clusters(clusters);

    // //merge compressed sets to each other and cluster
    // merge_compressed_sets(&C, &clusters);

    // //merge retained set to clusters
    // merge_retained_set(&R, &clusters);

    // //update centroids
    // update_centroids(&clusters);

    // free memory
    // for (i = 0; i < K; i++){
    //     free(clusters[i].centroid.coords);
    //     free(clusters[i].sum.coords);
    //     free(clusters[i].sum_squares.coords);
    // }
    free(clusters);

    // for (i = 0; i < R->number_of_points; i++){
    //     free(R->points[i].coords);
    // }
    free(R.points);
    // free(R);

    // for(i = 0; i < C->number_of_sets; i++){
    //     free(C->sets[i].points);
    // }
    // free(C->sets);
    // free(C);
    free(C.sets);
    fclose(input_file);
    return 0;
}