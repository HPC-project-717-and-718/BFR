#ifndef DEFINITIONS
#define DEFINITIONS

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <float.h> // needed for DBL_MAX

#define K 2                     // number of clusters
#define M 2                     // number of dimensions
#define S 100                   // sample size
#define T 3.5                   // threshold
#define LIMIT_S 10              // limit for rand generation of coordinates for initial centroids
#define MAX_SIZE_OF_BUFFER 1000 // max size of buffer in points
#define data_streamer FILE *    // data streamer type

#endif