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

ShiftK.out:$(OBJS) ../../src/libkomega.a
	$(F90) $(OBJS) -o $@ $(LIB) $(FFLAGS)

.F90.o:
	$(F90) -c $< $(FFLAGS) $(INCLUDE)

clean:
	rm -f *.o ShiftK.out *.mod

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
