MPI_PATH = /usr/local/bin

CC = $(MPI_PATH)/mpicc
FC = $(MPI_PATH)/mpif90

FC = ifort
CC = icc

OPT_FLAG = -g

CFLAGS  = $(OPT_FLAG)
FFLAGS  = $(OPT_FLAG)

OBJS = params.o FFTW.o Lai2D_1.o zgtsv.o xerbla.o
LIBS = -lfftw3 -lm

%.o : %.f
	$(FC) -c $< $(FFLAGS)

%.o : %.f90
	$(FC) -c $< $(FFLAGS)

%.o : %.F90
	$(FC) -c $< $(FFLAGS)


all : $(OBJS) poisson2D

poisson2D : $(OBJS) poisson2D.o
	$(FC) -o $@ $(LDFLAGS) poisson2D.o $(OBJS) $(LIBS)

poisson_c : $(OBJS) poisson_c.o
	$(CC) -o $@ $(LDFLAGS) poisson_c.o $(OBJS) $(LIBS)

poisson-multigrid : $(OBJS) poisson-multigrid.o
	$(FC) -o $@ $(LDFLAGS) poisson-multigrid.o $(OBJS) $(LIBS)

poisson2D.o : poisson2D.F90
Lai2D_1.o : Lai2D_1.F90

multigrid.o : multigrid.f90
	$(FC) -c $< $(FFLAGS)

jacobi.o : jacobi.f90

check : poisson2D
	@rm -f *.vti
	@poisson2D

clean :
	rm -rf $(OBJS) *.o *.vti *.mod poisson
