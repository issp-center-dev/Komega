MAINTAINERCLEANFILES    = Makefile.in
noinst_PROGRAMS = solve_rr.x solve_rc.x solve_cr.x solve_cc.x
solve_rr_x_SOURCES = solve_rr.F90
solve_rc_x_SOURCES = solve_rc.F90
solve_cr_x_SOURCES = solve_cr.F90
solve_cc_x_SOURCES = solve_cc.F90

AM_FCFLAGS = -I$(top_builddir) -I$(top_builddir)/src \
             -I$(top_srcdir) -I$(top_srcdir)/src
solve_rr_x_LDADD = -L $(top_builddir)/src -lkomega $(LAPACK_LIBS) $(BLAS_LIBS)
solve_rc_x_LDADD = -L $(top_builddir)/src -lkomega $(LAPACK_LIBS) $(BLAS_LIBS)
solve_cr_x_LDADD = -L $(top_builddir)/src -lkomega $(LAPACK_LIBS) $(BLAS_LIBS)
solve_cc_x_LDADD = -L $(top_builddir)/src -lkomega $(LAPACK_LIBS) $(BLAS_LIBS)

diagonalize.o : diagonalize.F90
make_ham.o : make_ham.F90 mathlib.o
mathlib.o : mathlib.F90
solve_cc.o : solve_cc.F90 mathlib.o
solve_cr.o : solve_cr.F90 mathlib.o
solve_rc.o : solve_rc.F90 mathlib.o
solve_rr.o : solve_rr.F90 mathlib.o
