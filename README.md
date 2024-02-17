# BFR

Github repository for the HPC course's project - Parallelization of the BFR algorithm

## Structure

```
.
├── README.md
├── data
│   ├── file1
│   └── synthetic
│       ├── synthetic_d10.txt
│       ├── synthetic_d2.txt
│       ├── synthetic_d2_4000points.txt
│       ├── synthetic_d2_4000points_2gaussians.txt
│       ├── synthetic_d3.txt
│       ├── synthetic_d3_6000points.txt
│       ├── synthetic_d4.txt
│       ├── synthetic_d5.txt
│       ├── synthetic_d6.txt
│       ├── synthetic_d7.txt
│       ├── synthetic_d8.txt
│       └── synthetic_d9.txt
├── lib
│   ├── bfr_structures
│   │   ├── bfr_structures.c
│   │   └── bfr_structures.h
│   ├── definitions.h
│   ├── kmeans
│   │   ├── LICENSE.md
│   │   ├── Makefile
│   │   ├── README.md
│   │   ├── example1.c
│   │   ├── example2.c
│   │   ├── kmeans.c
│   │   └── kmeans.h
│   └── kmeans_wrapper
│       ├── kmeans_wrapper.c
│       └── kmeans_wrapper.h
├── output
│   └── .gitkeep
├── src
│   ├── BFR_serial.c
│   └── bfr_serial
└── utils
    ├── 2dplot_d2_40.png
    ├── 2dplot_d2_4000.png
    ├── 3dplot_d3_60.png
    ├── 3dplot_d3_6000.png
    ├── create_synthetic_data.py
    ├── plot_data.py
    ├── pseudocode_BFR.md
    ├── testplot.png
    └── testplot2.png
```

The folders are organized as the following:

- ```data``` contains the datasets used by the BFR algorithm, with  ```synthetic``` being the folder containing data directly sampled from gaussians built with the data creation script
- ```output``` contains the algorithms' output file(s)
- ```src``` contains the algorithm's source code
- ```lib``` contains the source files' libraries
- ```utils``` contains useful files, notes and scripts used throughout the project's development cycle

## How to use

Compile the serial source code with:

```gcc src/BFR_serial.c lib/bfr_helper_functions/bfr_helper_functions.c lib/kmeans/kmeans.c lib/bfr_structures/bfr_structures.c  -o src/bfr_serial -lm```

Compile the parallel source code with:

```mpicc src/BFR_parallel.c lib/bfr_structures/bfr_structures.c lib/bfr_helper_functions/bfr_helper_functions.c lib/kmeans/kmeans.c -o src/BFR_parallel.o -lm -fopenmp```

