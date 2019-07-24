FORT = gfortran
CFLAGS = -c
LIBS = -llapack -lblas
#LIBS = lapack.f blas.f

all: delfan delneighbors
	./delfan 2 20 10 SAMPLE-2D-20N.dat
	./delneighbors 2 20 10 SAMPLE-2D-20N.dat

delfan: test_fan.f90 delaunayfan.o afl.o
	$(FORT) test_fan.f90 delaunayfan.o afl.o $(LIBS) -o delfan

delneighbors: test_neighbors.f90 delaunayneighbors.o afl.o
	$(FORT) test_neighbors.f90 delaunayneighbors.o afl.o $(LIBS) \
		 -o delneighbors

afl.o: AFL.f90
	$(FORT) $(CFLAGS) AFL.f90 -o afl.o

delaunayfan.o: delaunayfan.f90 afl.o
	$(FORT) $(CFLAGS) delaunayfan.f90 -o delaunayfan.o

delaunayneighbors.o: delaunayneighbors.f90 afl.o
	$(FORT) $(CFLAGS) delaunayneighbors.f90 -o delaunayneighbors.o

clean:
	rm -f *.o *.mod delfan delneighbors
