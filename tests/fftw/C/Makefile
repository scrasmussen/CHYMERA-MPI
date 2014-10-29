MPI_PATH = /usr/local/bin

CC = $(MPI_PATH)/mpicc
FC = $(MPI_PATH)/mpif90

FC = gfortran
CC = cc

OPT_FLAG = -debug -O0
OPT_FLAG = -g

CFLAGS  = $(OPT_FLAG) -I/usr/local/include
FFLAGS  = $(OPT_FLAG)
LDFLAGS = $(OPT_FLAG) -Wl,-no_pie -L/usr/local/lib

#CFLAGS +=  -I/usr/local/packages/fftw-3.3.4/api  -I/usr/local/packages/fftw-3.3.4/kernel  -I/usr/local/packages/fftw-3.3.4/rdft

OBJS = 
LIBS = -lfftw3 -lm

%.o : %.c
	$(CC) -c $< $(CFLAGS)

%.o : %.f
	$(FC) -c $< $(FFLAGS)

%.o : %.f90.
	$(FC) -c $< $(FFLAGS)

%.o : %.F90
	$(FC) -c $< $(FFLAGS)


all : $(OBJS) test_fftw

test_fftw : $(OBJS) test_fftw.o
	$(CC) -o $@ $(LDFLAGS) test_fftw.o $(OBJS) $(LIBS)

check : test_fftw
	@rm -f *.vti
	@test_fftw

clean :
	rm -rf $(OBJS) *.o test_fftw
	rm -rf *.mod *.i90 *.vti
