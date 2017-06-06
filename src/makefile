include ../make.sys

.SUFFIXES :
.SUFFIXES : .o .F90

OBJS = \
komega_cg_c.o \
komega_cg_r.o \
komega_cocg.o \
komega_bicg.o \
komega_math.o \
komega_vals.o

libshiftk:libkomega.a

libkomega.a:$(OBJS)
	ar cr libkomega.a $(OBJS)
#	cd shared; $(F90) -shared $(OBJS) $(FFLAGS) -fPIC -o libkomega.so

.F90.o:
	$(F90) -c $< $(FFLAGS)
#	cd shared;$(F90) -c ../$< $(FFLAGS) -fPIC

clean:
	rm -f *.o *.a *.mod

komega_cg_c.o:komega_math.o
komega_cg_c.o:komega_vals.o
komega_cg_r.o:komega_math.o
komega_cg_r.o:komega_vals.o
komega_cocg.o:komega_math.o
komega_cocg.o:komega_vals.o
komega_bicg.o:komega_math.o
komega_bicg.o:komega_vals.o
komega_math.o:komega_vals.o

