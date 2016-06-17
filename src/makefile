
F90=ifort
FFLAGS=-g -traceback -mkl -check all -warn all

.SUFFIXES :
.SUFFIXES : .o .f90

OBJS = \
mathlib.o \
shifted_krylov.o \
solve_rr.o

.f90.o:
	$(F90) -c $< $(FFLAGS)

solve_rr.x:$(OBJS)
	$(F90) $(OBJS) -o $@ $(FFLAGS)

clean:
	rm -f *.o *.x *.mod

mathlib.o:mathlib.f90
shifted_krylov.o:shifted_krylov.f90
shifted_krylov.o:mathlib.o
solve_rr.o:solve_rr.f90
solve_rr.o:shifted_krylov.o
solve_rr.o:mathlib.o

