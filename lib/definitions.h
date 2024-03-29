#ifndef DEFINITIONS
#define DEFINITIONS

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <float.h> // needed for DBL_MAX
#include <assert.h>
#include <string.h>

#define K 10                                    // number of clusters
#define M 2                                     // number of dimensions
#define S 100                                   // sample size
#define T 3.5                                   // threshold
#define LIMIT_S 10                              // limit for rand generation of coordinates for initial centroids
#define MAX_SIZE_OF_BUFFER 1000                 // max size of buffer in points
#define NORMAL_KMEANS_MAX_ITERATIONS 200        // max number of normal kmeans iterations
#define DEBUG 0                                 // turn on/off debug prints
#define DEBUG_TIME 1
#define data_streamer FILE *                    // data streamer type
#define BETA 4                                  // beta parameter for BFR algorithm
#define DBL_MAX_HC 10000.
#define K3 3

// #define KMEANS_THREADED                         // comment out for serial version

# define MASTER 0
# define NUMBER_OF_THREADS 4
# define DATA_BUFFER_SIZE 1000                  // equal to serial's MAX_SIZE_OF_BUFFER / NUMBER_OF_THREADS
# define MIN_DATA_BUFFER_SIZE_LAST_ROUND 50     // minimum number of points per node at the last round
# define UPPER_BOUND_ITERATIONS 200             // equal to serial's UPPER_BOUND_ITERATIONS

#endif