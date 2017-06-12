#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([komega], [1.1.0], [mkawamura@issp.u-tokyo.ac.jp])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_HEADERS([config.h])

echo ""
echo "  Check compiler"
echo ""
AC_PROG_F77
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
                                 [compile with MPI (parallelization) support.
                                  If none is found, MPI is not used. Default: no])],,
            [with_mpi=no])
if test x"$with_mpi" != xno; then
  AX_PROG_FC_MPI(,[FCFLAGS="$FCFLAGS -D__MPI"],[
    AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
    ])
else
  AC_PROG_FC
fi
echo ""
echo "  Check OpenMP"
echo ""
AC_ARG_WITH(openmp, [AS_HELP_STRING([--with-openmp],
                                 [compile with OpenMP (parallelization) support.
                                  If none is found, OpenMP is not used. Default: yes])],,
            [with_openmp=yes])
if test x"$with_openmp" != xno; then
  AX_OPENMP([FCFLAGS="$FCFLAGS $OPENMP_CFLAGS"],[
    AC_MSG_FAILURE([OpenMP requested, but couldn't use OpenMP.])
  ])
fi
echo ""
echo "  Check library generation"
echo ""
AM_PROG_AR
AC_PROG_LIBTOOL
echo ""
echo "  Check BLAS and LAPACK"
echo ""
AX_BLAS
AX_LAPACK
#
# Shread-safe
#
AC_ARG_ENABLE(threadsafe, [AS_HELP_STRING([--enable-threadsafe],
                                 [Thread-safe mode (experimental).])],,
                           [enable_threadsafe=no])
if test x"$enable_threadsafe" != xno; then
  FCFLAGS="$FCFLAGS -D__KOMEGA_THREAD"
fi
#
# zdotc and adotu
#
AC_ARG_ENABLE(zdot, [AS_HELP_STRING([--disable-zdot],
                           [Do not use ZDOTC and ZDOTU in BLAS.])],,
                           [enable_zdot=yes])
if test x"$enable_zdot" != xyes; then
  FCFLAGS="$FCFLAGS -D__NO_ZDOT"
fi
echo ""
echo "  Output files"
echo ""
AC_CONFIG_FILES([
Makefile
src/Makefile
test/Makefile
app/src/Makefile
])

AC_OUTPUT

echo ""
echo "  SUMMARY"
echo ""
echo "             FC = ${FC}"
echo "        FCFLAGS = ${FCFLAGS}"
echo "    LAPACK_LIBS = ${LAPACK_LIBS}"
echo "      BLAS_LIBS = ${BLAS_LIBS}"
echo ""
echo "  Configure done."
echo ""
