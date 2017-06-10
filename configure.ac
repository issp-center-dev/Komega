#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([komega], [1.1.0], [mkawamura@issp.u-tokyo.ac.jp])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_HEADERS([config.h])

# If --with-mpi=auto is used, try to find MPI, but use standard FC compiler if it is not found.
# If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
# If --with-mpi=no is used, use a standard FC compiler instead.
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
    [compile with MPI (parallelization) support. If none is found,
    MPI is not used. Default: no])
],,[with_mpi=no])
 AX_PROG_FC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
  use_mpi=no
  if test x"$with_mpi" = xyes; then
    AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
  else
    AC_MSG_WARN([No MPI compiler found, won't use MPI.])
  fi
])
AM_PROG_AR
AC_PROG_RANLIB
AX_BLAS
AX_LAPACK
AM_CONDITIONAL(USE_MPI, test x"$use_mpi" = xyes)

#AC_PROG_LIBTOOL

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
echo "           FC = ${FC}"
echo "      FCFLAGS = ${FCFLAGS}"
echo "  LAPACK_LIBS = ${LAPACK_LIBS}"
echo "    BLAS_LIBS = ${BLAS_LIBS}"
echo "          MPI = ${use_mpi}"
