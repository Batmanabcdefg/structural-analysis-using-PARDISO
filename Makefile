# (c) 2018, Jakub Mikula

# compiler
FC = ifort
MPI_FC = mpif90

FCFLAGS = -O3 -mkl -fopenmp
MPI_FCFLAGS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

# executable file
EXEC = kuf
MPI_EXEC = kuf_mpi

#..................................................................................

default: openmp

#shared memory architectures
openmp:
	$(MAKE) $(EXEC) $(FCFLAG)

#distributed memory arcitectures
mpi: 
	$(MAKE) $(MPI_EXEC) $(MPI_FCFLAGS)

#.................................................................................. 

$(EXEC) : main.f90
	$(FC) -o $(EXEC) main.f90 $(FCFLAGS) -fpp

$(MPI_EXEC) : main.f90
	$(MPI_FC) -o $(MPI_EXEC) main.f90 $(MPI_FCFLAGS) -cpp -DParallel

