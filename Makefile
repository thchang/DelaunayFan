FORT = gfortran
CFLAGS = -c
INC = -llapack -lblas

all: delfan

delfan: main.f90 delaunayfan.o afl.o vector.o
	$(FORT) main.f90 delaunayfan.o afl.o vector.o $(INC) -o delfan

afl.o: AFL.f90 vector.o
	$(FORT) $(CFLAGS) AFL.f90 -o afl.o

vector.o: Vector.f90
	$(FORT) $(CFLAGS) Vector.f90 -o vector.o

delaunayfan.o: delaunayfan.f90 afl.o
	$(FORT) $(CFLAGS) delaunayfan.f90 -o delaunayfan.o

clean:
	rm *.o *.mod vtdel
