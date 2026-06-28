MAINTAINERCLEANFILES    = Makefile.in
noinst_PROGRAMS = solve_rr.x solve_rc.x solve_cr.x solve_cc.x
solve_rr_x_SOURCES = solve_rr.F90
solve_rc_x_SOURCES = solve_rc.F90
solve_cr_x_SOURCES = solve_cr.F90
solve_cc_x_SOURCES = solve_cc.F90

AM_FCFLAGS = -I$(top_builddir) -I$(top_builddir)/src 
solve_rr_x_LDADD = -L$(top_builddir)/src -lkomega $(LAPACK_LIBS) $(BLAS_LIBS)
solve_rc_x_LDADD = -L$(top_builddir)/src -lkomega $(LAPACK_LIBS) $(BLAS_LIBS)
solve_cr_x_LDADD = -L$(top_builddir)/src -lkomega $(LAPACK_LIBS) $(BLAS_LIBS)
solve_cc_x_LDADD = -L$(top_builddir)/src -lkomega $(LAPACK_LIBS) $(BLAS_LIBS)

diagonalize.o : diagonalize.F90
make_ham.o : make_ham.F90 mathlib.o
mathlib.o : mathlib.F90
solve_cc.o : solve_cc.F90 mathlib.o
solve_cr.o : solve_cr.F90 mathlib.o
solve_rc.o : solve_rc.F90 mathlib.o
solve_rr.o : solve_rr.F90 mathlib.o

# Helper sources compiled via the explicit rules below (used by the drivers),
# and the regression tests run by "make check" (see run_tests.sh). All must be
# distributed so the test directory builds and "make check" works from a
# release tarball.
EXTRA_DIST = mathlib.F90 make_ham.F90 diagonalize.F90 \
             run_tests.sh test_lzconv.F90 \
             regression_real.in regression_complex.in \
             check_configure_zdot.sh

check-local: all-am
	srcdir='$(srcdir)' top_builddir='$(top_builddir)' \
	FC='$(FC)' FCFLAGS='$(FCFLAGS)' \
	LAPACK_LIBS='$(LAPACK_LIBS)' BLAS_LIBS='$(BLAS_LIBS)' \
	LIBTOOL='$(top_builddir)/libtool' \
	$(SHELL) $(srcdir)/run_tests.sh
