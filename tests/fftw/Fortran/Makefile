FC = ifort
CC = icc

OPT_FLAG = -debug -O0

CFLAGS  = $(OPT_FLAG)
FFLAGS  = $(OPT_FLAG)
LDFLAGS = $(OPT_FLAG) -Wl,-no_pie -L/usr/local/lib

OBJS = params.o FFTW.o DFT.o
LIBS = -lfftw3 -lm

%.o : %.f
	$(FC) -c $< $(FFLAGS)

%.o : %.f90
	$(FC) -c $< $(FFLAGS)

%.o : %.F90
	$(FC) -c $< $(FFLAGS)

all : $(OBJS) test_DFT

test_DFT.o : test_DFT.F90
test_DFT : $(OBJS) test_DFT.o
	$(FC) -o $@ $(LDFLAGS) test_DFT.o $(OBJS) $(LIBS)

check : test_DFT
	@test_DFT

clean :
	rm -rf $(OBJS) *.o test_DFT
	rm -rf *.mod
