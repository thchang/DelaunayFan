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

### Prerequisites

Before installing the DelaunayFan subroutine, you will need BLAS and LAPACK,
which are available on NIST - Guide to Available Math Software (GAMS) and
Netlib.
Alternatively, use any package manager to install
```
libblas
```
and
```
liblapack
```

### Contents

There are 9 files included:
 - delaunayfan.f90 contains the main subroutine DelaunayFan.
 - AFL.f90 contains the source code for the data structure described in Chang et. al.
 - Vector.f90 contains the source code for a dynamic multidimensional array, similar to the C++ vector objects.
 - main.f90 contains driver code generating a command line program for DelaunayFan
 - Makefile for building the code.
 - SAMPLE-2D-20N.dat contains a pseudo-randomly generated input set in 2-dimensions containing 20 vertices/input points.
 - LICENSE contains a description of the MIT license.
 - README.md is this file.
 - delvor-ieeesec.pdf contains a detailed description of the algorithm.

### Compiling and Testing

Use the included Makefile to compile the main program and all the subroutines:
```
make -B
```
Then, try running the sample executable with the provided input file:
```
./delfan 2 20 [any integer from 1-20] SAMPLE-2D-20N.dat
```

## License

This project is licensed under the MIT License. See the LICENSE file.

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
