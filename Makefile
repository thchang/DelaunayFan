FORT = gfortran
CFLAGS = -c
INC = -llapack -lblas

all: delfan

delfan: main.f90 delaunayfan.o afl.o
	$(FORT) main.f90 delaunayfan.o afl.o $(INC) -o delfan

afl.o: AFL.f90
	$(FORT) $(CFLAGS) AFL.f90 -o afl.o

delaunayfan.o: delaunayfan.f90 afl.o
	$(FORT) $(CFLAGS) delaunayfan.f90 -o delaunayfan.o

clean:
	rm *.o *.mod delfan
