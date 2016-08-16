include ../make.sys

.SUFFIXES :
.SUFFIXES : .o .f90

LIB = -L ../src -lshiftk mathlib.o
INCLUDE = -I ../src

all:solve_rr.x solve_cr.x solve_rc.x solve_cc.x

solve_cr.x:solve_cr.o ../src/libshiftk.a
	$(F90) $< -o $@ $(FFLAGS) $(LIB)

solve_rr.x:solve_rr.o ../src/libshiftk.a
	$(F90) $< -o $@ $(FFLAGS) $(LIB)

solve_rc.x:solve_rc.o ../src/libshiftk.a
	$(F90) $< -o $@ $(FFLAGS) $(LIB)

solve_cc.x:solve_cc.o ../src/libshiftk.a
	$(F90) $< -o $@ $(FFLAGS) $(LIB)

.f90.o:
	$(F90) -c $< $(FFLAGS) $(INCLUDE)

clean:
	rm -f *.o *.x *.mod

solve_cr.o:mathlib.o
solve_cr.o:../src/libshiftk.a
solve_rr.o:mathlib.o
solve_rr.o:../src/libshiftk.a
solve_rc.o:mathlib.o
solve_rc.o:../src/libshiftk.a
solve_cc.o:mathlib.o
solve_cc.o:../src/libshiftk.a

