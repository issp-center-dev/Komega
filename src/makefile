include ../make.sys

.SUFFIXES :
.SUFFIXES : .o .f90

OBJS = \
shifted_krylov.o \
shifted_krylov_math.o \
shifted_krylov_vals.o

libshiftk.a:$(OBJS)
	ar cr $@ $(OBJS)

.f90.o:
	$(F90) -c $< $(FFLAGS)

clean:
	rm -f *.o *.a *.mod

shifted_krylov.o:shifted_krylov_math.o
shifted_krylov.o:shifted_krylov_vals.o

