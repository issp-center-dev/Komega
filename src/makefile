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

libshiftk:libshiftk.a shared/libshiftk.so

libshiftk.a:$(OBJS)
	ar cr $@ $(OBJS)

shared/libshiftk.so:$(OBJS)
	mkdir shared
	$(F90) -shared $(OBJS) $(FFLAGS) -fPIC -o $@

.f90.o:
	$(F90) -c $< $(FFLAGS) -fPIC

clean:
	rm -f *.o *.a *.mod
	rm -r -f shared

shifted_cg_c.o:shifted_krylov_math.o
shifted_cg_c.o:shifted_krylov_vals.o
shifted_cg_r.o:shifted_krylov_math.o
shifted_cg_r.o:shifted_krylov_vals.o
shifted_cocg.o:shifted_krylov_math.o
shifted_cocg.o:shifted_krylov_vals.o
shifted_bicg.o:shifted_krylov_math.o
shifted_bicg.o:shifted_krylov_vals.o
shifted_krylov_math.o:shifted_krylov_vals.o

