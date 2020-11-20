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
 - delaunayneighbors.f90 contains a slight modification to DELAUNAYFAN,
   DELAUNAYNEIGHBORS, which returns the indices of the Delaunay neighbors,
   but not the simplices
 - AFL.f90 contains the source code for the face list data structure described in Chang et. al, as well as a subset of the Vector class (see thchang/VectorClass)
 - blas.f contains a minimal copy of BLAS, containing only the procedures required by delaunayfan.f90
 - lapack.f contains a minimal copy of LAPACK, containing only the procedures required by delaunayfan.f90
 - test\_fan.f90 contains driver code generating a command line program for
   DELAUNAYFAN
 - test\_neighbors.f90 contains driver code generating a command line program
   for DELAUNAYNEIGHBORS
 - Makefile for building the code using the gfortran compiler
 - SAMPLE-2D-20N.dat contains a pseudo-randomly generated input set in
   2-dimensions containing 20 vertices/input points
 - SAMPLE-4D-43N-DEGEN.dat contains a degenerate dataset from a real world
   problem in 4D with 43 data points. DELAUNAYFAN and DELAUNAYNEIGHBORS
   throw errors (IERR=40) when run on these datasets.
 - delvor-ieeesec.pdf contains a detailed description of the algorithm

### Compiling and Testing

Use the included Makefile to compile the main program and all the subroutines (requires gfortran):
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

```
@inproceedings{chang2018computing,
   author={T. H. Chang and L. T. Watson and T. C. H. Lux and S. Raghvendra and B. Li and L. Xu and A. R. Butt and K. W. Cameron and Y. Hong},
   booktitle={Proceedings of SoutheastCon 2018},
   title={Computing the Umbrella Neighbourhood of a Vertex in the {D}elaunay Triangulation and a Single {V}oronoi Cell in Arbitrary Dimension},
   year={2018},
   pages={1-8},
   publisher = {IEEE},
   location = {St. Petersburg, FL, USA},
   doi={10.1109/SECON.2018.8479003},
   month={April}
}
```

## Related Work

If you found this work interesting, you may also consider
https://github.com/vtopt/DelaunaySparse

