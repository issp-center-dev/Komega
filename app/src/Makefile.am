MAINTAINERCLEANFILES    = Makefile.in

bin_PROGRAMS = ShiftK.out

ShiftK_out_SOURCES = dyn_mod.F90 ham_prod.F90 lobpcg_mod.F90 \
                     shiftk.F90 shiftk_io.F90 shiftk_vals.F90

AM_FCFLAGS = -I$(top_builddir) -I$(top_builddir)/src \
             -I$(top_srcdir) -I$(top_srcdir)/src
if USE_MPI
  AM_FCFLAGS += -D __MPI
endif
ShiftK_out_LDADD = -L $(top_builddir)/src -lkomega $(LAPACK_LIBS) $(BLAS_LIBS)

dyn_mod.o : dyn_mod.F90 ham_prod.o lobpcg_mod.o shiftk_vals.o shiftk_io.o
ham_prod.o : ham_prod.F90 shiftk_vals.o
lobpcg_mod.o : lobpcg_mod.F90 ham_prod.o shiftk_vals.o
shiftk.o : shiftk.F90 dyn_mod.o ham_prod.o lobpcg_mod.o shiftk_vals.o shiftk_io.o
shiftk_io.o : shiftk_io.F90 lobpcg_mod.o ham_prod.o ham_prod.o shiftk_vals.o
shiftk_vals.o : shiftk_vals.F90
