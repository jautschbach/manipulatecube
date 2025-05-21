#
# This program may need lapack and blas routines
#
# For MPI support, we obviously also need the MPI libs and include files
#
# Note: we define two compiler variables: $FC for fixed format F77 or F90 and $FC90 for
# free format F90


### Fortran compiler
FC     = gfortran
FC90   = gfortran

LINKER = gfortran

### Disable/Enable OpenMP support
OMP = 

### Compilation/Linking options
FOPT = -O2

# Debug flags (used for 'make debug')
# gfortran:
FDBG = -g -fimplicit-none -fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow -Wall
# Intel fortran:
# FDBG = -g -traceback -C -fpe0

### Extra libraries
FLIB=-llapack   # use this line if you have Lapack libraries installed
#FLIB =         # use this line if you're using source code for DLASRT

### Include dir(s)
INC =

FFLAGS = $(OMP) $(INC)

.SUFFIXES: .f .F .F90 .f90

BIN = manipulatecube

OBJ =      # use this line if you have Lapack libraries installed
#OBJ = xerbla.o lsame.o dlasrt.o # use this when using source code for DLASRT

OBJ90 = types.o 


.F.o:
	$(FC) $(FOPT) $(FFLAGS) -c $< -o $@
.F90.o: 
	$(FC90) $(FOPT) $(FFLAGS) -c $< -o $@
.f90.o: 
	$(FC90) $(FOPT) $(FFLAGS) -c $< -o $@

all: $(BIN)

manipulatecube: $(OBJ) $(OBJ90) manipulatecube.o
	$(LINKER) -o $(BIN).exe manipulatecube.o $(OBJ) $(OBJ90) $(FFLAGS) $(FLIB)

clean:
	rm -f *.o *.mod

realclean:
	rm -f $(BIN).exe *.o *.mod

debug: FFLAGS += $(FDBG)
debug: all

