include ../../make.sys

.SUFFIXES :
.SUFFIXES : .o .F90

LIB = -L ../../src/ -lkomega
MPILIB = -L ../../../src/mpi/ -lpkomega
INCLUDE = -I ../../src/
MPIINCLUDE = -I ../../../src/mpi/

OBJS= \
dyn_mod.o \
ham_prod.o \
lobpcg_mod.o \
shiftk.o \
shiftk_io.o \
shiftk_vals.o \

ifdef MPIF90
ShiftK.out:$(OBJS) ../../src/libkomega.a
	$(F90) $(OBJS) -o $@ $(LIB) $(FFLAGS)
	cd mpi; $(MPIF90) $(OBJS) -o $@ $(MPILIB) $(FFLAGS) -D MPI
else
ShiftK.out:$(OBJS) ../../src/libkomega.a
	$(F90) $(OBJS) -o $@ $(LIB) $(FFLAGS)
endif

ifdef MPIF90
.F90.o:
	$(F90) -c $< $(FFLAGS) $(INCLUDE)
	cd mpi; $(MPIF90) -c ../$< $(FFLAGS) $(MPIINCLUDE) -D MPI
else
.F90.o:
	$(F90) -c $< $(FFLAGS) $(INCLUDE)
endif

clean:
	rm -f *.o ShiftK.out *.mod
	cd mpi; rm -f *.o ShiftK.out *.mod

shiftk.o:shiftk_io.o
shiftk.o:dyn_mod.o
shiftk.o:ham_prod.o
shiftk.o:lobpcg_mod.o
ham_prod.o:shiftk_vals.o
lobpcg_mod.o:shiftk_vals.o
lobpcg_mod.o:ham_prod.o
shiftk_io.o:shiftk_vals.o
shiftk_io.o:ham_prod.o
shiftk_io.o:lobpcg_mod.o
dyn_mod.o:../../src/libkomega.a
dyn_mod.o:shiftk_io.o
dyn_mod.o:shiftk_vals.o
dyn_mod.o:ham_prod.o
dyn_mod.o:lobpcg_mod.o
