include ../make.sys

.SUFFIXES :
.SUFFIXES : .o .f90

LIB = -L ../src -lshiftk
INCLUDE = -I ../src

solve_rr.x:solve_rr.o 
	$(F90) $< -o $@ $(FFLAGS) $(LIB)

.f90.o:
	$(F90) -c $< $(FFLAGS) $(INCLUDE)

clean:
	rm -f *.o *.x *.mod

solve_rr.o:mathlib.o

