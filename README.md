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
│       ├── synthetic_d3.txt
│       ├── synthetic_d3_6000points.txt
│       ├── synthetic_d4.txt
│       ├── synthetic_d5.txt
│       ├── synthetic_d6.txt
│       ├── synthetic_d7.txt
│       ├── synthetic_d8.txt
│       └── synthetic_d9.txt
├── output
├── src
│   └── BFR_serial.c
└── utils
    ├── 2dplot_d2_40.png
    ├── 2dplot_d2_4000.png
    ├── 3dplot_d3_60.png
    ├── 3dplot_d3_6000.png
    ├── create_synthetic_data.py
    ├── plot_data.py
    └── pseudocode_BFR.md
```

The folders are organized as the following:

- ```data``` contains the datasets used by the BFR algorithm, with  ```synthetic``` being the folder containing data directly sampled from gaussians built with the data creation script
- ```output``` contains the algorithms' output file(s)
- ```src``` contains the algorithm's source code
- ```utils``` contains useful files, notes and scripts used throughout the project's development cycle
