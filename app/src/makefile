include ../../make.sys

.SUFFIXES :
.SUFFIXES : .o .f90

LIB = -L ../../src/ -lshiftk
INCLUDE = -I ../../src/

all:shiftk.x

shiftk.x:shiftk.o ../../src/libshiftk.a
	$(F90) $< -o $@ $(FFLAGS) $(LIB)

.f90.o:
	$(F90) -c $< $(FFLAGS) $(INCLUDE)

clean:
	rm -f *.o *.x *.mod

shiftk.o:mathlib.o
shiftk.o:../../src/libshiftk.a

