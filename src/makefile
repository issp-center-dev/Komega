include ../make.sys

.SUFFIXES :
.SUFFIXES : .o .f90

OBJS = \
shifted_cg_c.o \
shifted_cg_r.o \
shifted_cocg.o \
shifted_bicg.o \
shifted_krylov_math.o \
shifted_krylov_vals.o

libshiftk.a:$(OBJS)
	ar cr $@ $(OBJS)

.f90.o:
	$(F90) -c $< $(FFLAGS)

clean:
	rm -f *.o *.a *.mod

shifted_cg_c.o:shifted_krylov_math.o
shifted_cg_c.o:shifted_krylov_vals.o
shifted_cg_r.o:shifted_krylov_math.o
shifted_cg_r.o:shifted_krylov_vals.o
shifted_cocg.o:shifted_krylov_math.o
shifted_cocg.o:shifted_krylov_vals.o
shifted_bicg.o:shifted_krylov_math.o
shifted_bicg.o:shifted_krylov_vals.o

