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

#define K 10                              // number of clusters
#define M 2                              // number of dimensions
#define DIMENSION 2                      // number of dimensions
#define S 100                            // sample size
#define T 3.5                            // threshold
#define LIMIT_S 10                       // limit for rand generation of coordinates for initial centroids
#define MAX_SIZE_OF_BUFFER 1000          // max size of buffer in points
#define NORMAL_KMEANS_MAX_ITERATIONS 200 // max number of normal kmeans iterations
#define DEBUG 1                          // turn on/off debug prints
#define data_streamer FILE *             // data streamer type
#define BETA 4                           // beta parameter for BFR algorithm

#define KMEANS_THREADED                  // comment out for serial version

# define MASTER 0
# define NUMBER_OF_THREADS 4
# define DATA_BUFFER_SIZE 250 // equal to serial's MAX_SIZE_OF_BUFFER / NUMBER_OF_THREADS

#endif