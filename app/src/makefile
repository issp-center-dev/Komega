include ../../make.sys

.SUFFIXES :
.SUFFIXES : .o .f90

LIB = -L ../../src/ -lshiftk
INCLUDE = -I ../../src/

OBJS= \
dyn_mod.o \
ham_prod.o \
lobpcg_mod.o \
shiftk.o \
shiftk_io.o \
shiftk_vals.o \

all:ShiftK.out

ShiftK.out:$(OBJS) ../../src/libshiftk.a
	$(F90) $(OBJS) -o $@ $(LIB) $(FFLAGS)

.f90.o:
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
dyn_mod.o:../../src/libshiftk.a
dyn_mod.o:shiftk_io.o
dyn_mod.o:shiftk_vals.o
dyn_mod.o:ham_prod.o
