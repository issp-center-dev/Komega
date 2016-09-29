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

libshiftk:libshiftk.a

ifdef MPIF90
libshiftk.a:$(OBJS)
	ar cr libshiftk.a $(OBJS)
	cd shared; $(F90) -shared $(OBJS) $(FFLAGS) -fPIC -o libshiftk.so
	cd mpi; ar cr libpshiftk.a $(OBJS)
	cd shared_mpi; $(F90) -shared $(OBJS) $(FFLAGS) -fPIC -o libpshiftk.so -D MPI
else
libshiftk.a:$(OBJS)
	ar cr libshiftk.a $(OBJS)
	cd shared; $(F90) -shared $(OBJS) $(FFLAGS) -fPIC -o libshiftk.so
endif

ifdef MPIF90
.f90.o:
	$(F90) -c $< $(FFLAGS)
	cd shared; $(F90) -c ../$< $(FFLAGS) -fPIC
	cd mpi; $(MPIF90) -c ../$< $(FFLAGS) -D MPI
	cd shared_mpi; $(MPIF90) -c ../$< $(FFLAGS) -fPIC -D MPI
else
.f90.o:
	$(F90) -c $< $(FFLAGS)
	cd shared;$(F90) -c ../$< $(FFLAGS) -fPIC
endif

clean:
	rm -f *.o *.a *.mod
	cd shared; rm -f *.o *.so *.mod
	cd mpi; rm -f *.o *.a *.mod
	cd shared_mpi; rm -f *.o *.so *.mod

shifted_cg_c.o:shifted_krylov_math.o
shifted_cg_c.o:shifted_krylov_vals.o
shifted_cg_r.o:shifted_krylov_math.o
shifted_cg_r.o:shifted_krylov_vals.o
shifted_cocg.o:shifted_krylov_math.o
shifted_cocg.o:shifted_krylov_vals.o
shifted_bicg.o:shifted_krylov_math.o
shifted_bicg.o:shifted_krylov_vals.o
shifted_krylov_math.o:shifted_krylov_vals.o

