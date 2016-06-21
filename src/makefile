include ../make.sys

.SUFFIXES :
.SUFFIXES : .o .f90

OBJS = \
mathlib.o \
shifted_krylov.o

libshiftk.a:$(OBJS)
	ar cr $@ $(OBJS)

.f90.o:
	$(F90) -c $< $(FFLAGS)

clean:
	rm -f *.o *.a *.mod

mathlib.o:mathlib.f90
shifted_krylov.o:shifted_krylov.f90
shifted_krylov.o:mathlib.o

