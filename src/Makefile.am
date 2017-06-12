#
# Komega library
#
MAINTAINERCLEANFILES = Makefile.in

pkginclude_HEADERS = komega.h 
lib_LTLIBRARIES = libkomega.la

libkomega_la_SOURCES = komega_bicg.F90 komega_cg_r.F90 komega_cg_c.F90 komega_cocg.F90 \
                       komega_math.F90 komega_vals.F90

komega_bicg.lo : komega_bicg.F90 komega_math.lo komega_vals.lo komega_vals.lo komega_vals.lo 
komega_cg_c.lo : komega_cg_c.F90 komega_math.lo komega_vals.lo komega_vals.lo komega_vals.lo 
komega_cg_r.lo : komega_cg_r.F90 komega_math.lo komega_vals.lo komega_vals.lo komega_vals.lo 
komega_cocg.lo : komega_cocg.F90 komega_math.lo komega_vals.lo komega_vals.lo komega_vals.lo 
komega_math.lo : komega_math.F90 komega_vals.lo 
komega_vals.lo : komega_vals.F90 
