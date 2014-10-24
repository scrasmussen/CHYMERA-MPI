MPI_PATH = /usr/local/bin

CC = $(MPI_PATH)/mpicc
FC = $(MPI_PATH)/mpif90

FC = ifort
CC = icc

OPT_FLAG = -g
OPT_FLAG = -debug -O0

CFLAGS  = $(OPT_FLAG)
FFLAGS  = $(OPT_FLAG)
LDFLAGS = $(OPT_FLAG) -Wl,-no_pie -L/usr/local/lib

OBJS = params.o FFTW.o DFT.o zgtsv.o xerbla.o
LIBS = -lfftw3 -lm

%.o : %.f
	$(FC) -c $< $(FFLAGS)

%.o : %.f90
	$(FC) -c $< $(FFLAGS)

%.o : %.F90
	$(FC) -c $< $(FFLAGS)


all : $(OBJS) test_DFT

test_DFT : $(OBJS) test_DFT.o
	$(FC) -o $@ $(LDFLAGS) test_DFT.o $(OBJS) $(LIBS)

poisson_c : $(OBJS) poisson_c.o
	$(CC) -o $@ $(LDFLAGS) poisson_c.o $(OBJS) $(LIBS)

poisson-multigrid : $(OBJS) poisson-multigrid.o
	$(FC) -o $@ $(LDFLAGS) poisson-multigrid.o $(OBJS) $(LIBS)

test_DFT.o : test_DFT.F90
Lai2D_1.o : Lai2D_1.F90

multigrid.o : multigrid.f90
	$(FC) -c $< $(FFLAGS)

jacobi.o : jacobi.f90

check : test_DFT
	@rm -f *.vti
	@test_DFT

clean :
	rm -rf $(OBJS) *.o test_DFT
	rm -rf *.mod *.i90 *.vti
