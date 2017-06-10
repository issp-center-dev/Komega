#
# Komega library
#
MAINTAINERCLEANFILES    = Makefile.in

pkginclude_HEADERS      = komega_bicg.h komega_cg_c.h komega_cg_r.h komega_cocg.h 
#lib_LTLIBRARIES = libkomega.la
lib_LIBRARIES = libkomega.a
#libkomega_la_SOURCES = komega_bicg.F90 komega_cg_r.F90 komega_cg_c.F90 komega_cocg.F90
libkomega_a_SOURCES = komega_bicg.F90 komega_cg_r.F90 komega_cg_c.F90 komega_cocg.F90 \
                      komega_math.F90 komega_vals.F90
if USE_MPI
  AM_FCFLAGS = -D __MPI
endif

komega_bicg.o : komega_bicg.F90 komega_math.o komega_vals.o komega_vals.o komega_vals.o 
komega_cg_c.o : komega_cg_c.F90 komega_math.o komega_vals.o komega_vals.o komega_vals.o 
komega_cg_r.o : komega_cg_r.F90 komega_math.o komega_vals.o komega_vals.o komega_vals.o 
komega_cocg.o : komega_cocg.F90 komega_math.o komega_vals.o komega_vals.o komega_vals.o 
komega_math.o : komega_math.F90 komega_vals.o 
komega_vals.o : komega_vals.F90 
