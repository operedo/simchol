CC	= mpicc
CPP	= mpic++
OMP	= -openmp 
CFLAGS	= -O3 $(OMP)
LIBS  = -lm

all: simchol 

simchol.o: simchol.cpp
	$(CPP) $(CFLAGS) -c -o simchol.o simchol.cpp

simchol: simchol.o
	$(CPP) $(CFLAGS) -o simchol simchol.o $(LIBS)

mrproper: clean
	rm -f job_*  salida* 

clean:
	rm -rf *.o simchol 
