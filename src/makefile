include ../make.sys

.SUFFIXES :
.SUFFIXES : .o .f90

OBJS = \
komega_cg_c.o \
komega_cg_r.o \
komega_cocg.o \
komega_bicg.o \
komega_math.o \
komega_vals.o

libshiftk:libkomega.a

ifdef MPIF90
libkomega.a:$(OBJS)
	ar cr libkomega.a $(OBJS)
	cd shared; $(F90) -shared $(OBJS) $(FFLAGS) -fPIC -o libkomega.so
	cd mpi; ar cr libpkomega.a $(OBJS)
	cd shared_mpi; $(F90) -shared $(OBJS) $(FFLAGS) -fPIC -o libpkomega.so -D MPI
else
libkomega.a:$(OBJS)
	ar cr libkomega.a $(OBJS)
	cd shared; $(F90) -shared $(OBJS) $(FFLAGS) -fPIC -o libkomega.so
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

komega_cg_c.o:komega_math.o
komega_cg_c.o:komega_vals.o
komega_cg_r.o:komega_math.o
komega_cg_r.o:komega_vals.o
komega_cocg.o:komega_math.o
komega_cocg.o:komega_vals.o
komega_bicg.o:komega_math.o
komega_bicg.o:komega_vals.o
komega_math.o:komega_vals.o

