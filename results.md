# Platform: Darwin at LANL
 * RAM: 125GB	
 * CPU: Intel broadwell	E5-2695_v4,	2.10GHz
 * Interconnect: Infiniband, edr 2, 1 Gb

## File:
 * Argonne_L360_HACC001/STEP499/m000.full.mpicosmo.499
 * 3D Split: 8, 8, 4, physical coordinates: (0, 0, 0) -> (256, 256, 256) 
 * Ranks in file: 256, Partitions: 64
 * Total # particles: 1,073,726,359
 * Size: 46 G (from dh -h on folder STEP 499)

## Job size:
 64 MPI ranks (8 nodes with 8 mpi ranks per node)

## Outputing decompressed:
 * BigCrunch: Total run time 251.184 s, Size: 38 GB (one file)
 * SZ:        Total run time 222.912 s, Size: 38 GB (one file)


## Comptression metrics
| Compressor_field |  absolute_error |  relative_error |  mse        |  Max Com Throughput MB/s|  Max DeCom Throughput MB/s|  Min Com Throughput MB/s |  Min DeCom Throughput MB/s|  Compression Ratio | 
|------------------|-----------------|-----------------|-------------|-------------------------|---------------------------|--------------------------|---------------------------|--------------------| 
| BigCrunch_x      |  0.125          |  0.000976562    |  0.00287934 |  14.786                 |  117.611                  |  9.68785                 |  52.416                   |  13.553            |
| SZ_x             |  0.255783       |  0.000999547    |  0.00710903 |  50.2962                |  112.736                  |  37.4287                 |  53.1926                  |  16.5391           |  
| BigCrunch_y      |  0.125          |  0.000976562    |  0.00301344 |  14.4318                |  103.803                  |  10.2699                 |  46.2856                  |  13.0411           | 
| SZ_y             |  0.255814       |  0.000999548    |  0.0075367  |  49.9475                |  109.709                  |  33.9866                 |  49.4052                  |  15.3533           | 
| BigCrunch_z      |  0.125          |  0.000976562    |  0.00301163 |  9.78527                |  85.9057                  |  8.52105                 |  49.3865                  |  11.1931           | 
| SZ_z             |  0.255692       |  0.000998957    |  0.00705831 |  47.397                 |  88.3291                  |  36.0938                 |  73.9745                  |  13.2241           | 
| BigCrunch_vx     |  2              |  0.000976562    |  0.0162471  |  11.2349                |  50.1134                  |  8.82869                 |  28.706                   |  2.70519           | 
| SZ_vx            |  3.52686        |  0.000999044    |  0.0313098  |  42.1682                |  28.6738                  |  34.888                  |  25.829                   |  2.89944           | 
| BigCrunch_vy     |  2              |  0.000976562    |  0.0174506  |  11.2901                |  48.2419                  |  8.78459                 |  28.162                   |  2.72096           | 
| SZ_vy            |  3.69727        |  0.000999049    |  0.03365    |  42.0925                |  28.1931                  |  34.57                   |  25.7233                  |  2.91714           | 
| BigCrunch_vz     |  2              |  0.000976562    |  0.0188472  |  11.0076                |  46.9743                  |  8.70787                 |  30.515                   |  2.72669           | 
| SZ_vz            |  3.76123        |  0.000999046    |  0.0363292  |  42.2357                |  28.4567                  |  32.8614                 |  25.921                   |  2.92593           | 
