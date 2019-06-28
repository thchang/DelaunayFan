# Delaunay Fan

Compute the umbrella neighbourhood of a single vertex in the Delaunay
triangulation or (equivalently) a single Voronoi cell, without
computing the entire Delaunay triangulation.
Designed for usage in arbitrary dimension.
This research code is written in Fortran 2008 and is based on the algorithm
described in:

T. H. Chang, L. T. Watson, T. C.H. Lux, S. Raghvendra, B. Li, L. Xu,
A. R. Butt, K. W. Cameron, and Y. Hong. Computing the umbrella
neighbourhood of a vertex in the Delaunay triangulation and a single
Voronoi cell in arbitrary dimension. In Proceedings of IEEE Southeast
Conference 2018 (SEC 2018). IEEE, St. Petersburg, FL, 2018.

Because some degeneracy issues could not be resolved, this code is not
robust for degenerate or near-degenerate inputs.

## Getting Started

### Dependencies

DELAUNAYFAN depends on several procedures from BLAS and LAPACK.
For optimal performance, link a system installation using
```
-lblas
```
and
```
-llapack
```
Alternatively, minimal copies of BLAS and LAPACK are provided in the files
blas.f and lapack.f.

### Contents

 - delaunayfan.f90 contains the main subroutine DELAUNAYFAN
 - AFL.f90 contains the source code for the face list data structure described in Chang et. al, as well as a subset of the Vector class (see thchang/VectorClass)
 - blas.f contains a minimal copy of BLAS, containing only the procedures required by delaunayfan.f90
 - lapack.f contains a minimal copy of LAPACK, containing only the procedures required by delaunayfan.f90
 - main.f90 contains driver code generating a command line program for DelaunayFan
 - Makefile for building the code
 - SAMPLE-2D-20N.dat contains a pseudo-randomly generated input set in 2-dimensions containing 20 vertices/input points
 - delvor-ieeesec.pdf contains a detailed description of the algorithm

### Compiling and Testing

Use the included Makefile to compile the main program and all the subroutines:
```
make -B
```
Then, try running the sample executable with the provided input file:
```
./delfan 2 20 [any integer from 1-20] SAMPLE-2D-20N.dat
```

## Authors

* **Tyler H. Chang** - *Primary author*

For other authors, see the reference.

## Reference

To cite this work, use:

T. H. Chang, L. T. Watson, T. C.H. Lux, S. Raghvendra, B. Li, L. Xu,
A. R. Butt, K. W. Cameron, and Y. Hong. Computing the umbrella
neighbourhood of a vertex in the Delaunay triangulation and a single
Voronoi cell in arbitrary dimension. In Proceedings of IEEE Southeast
Conference 2018 (SEC 2018). IEEE, St. Petersburg, FL, 2018.
