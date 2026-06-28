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
      AC_MSG_WARN([OpenMP requested, but couldn't use OpenMP.])
  ])
fi
echo ""
echo "  Check fortran module"
echo ""
# Find out how to add to the Fortran include path.
AX_F90_MODULE_FLAG
# The module flag may be "-I " or so, and some make implementations
# swallow trailing space in variable settings, so we employ a little
# trick.
AC_SUBST([empty], [""])
AC_SUBST([FC_MODINC], [$ax_cv_f90_modflag'${empty}'])
# Find out the file name extensions of Fortran module files.
AX_F90_MODULE_EXTENSION
if test "x$ax_f90_modext" != xunknown; then
  FC_MODEXT=$ax_cv_f90_modext
else
  FC_MODEXT=$ax_cv_f90_modext
fi
AC_SUBST([FC_MODEXT])
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
# zdotc and zdotu
#
# The reference (f2c/g77) BLAS returns a COMPLEX result through a hidden first
# pointer argument, whereas komega_math.F90 declares zdotc/zdotu as ordinary
# register-return COMPLEX(8) FUNCTIONs. Calling such a BLAS through that
# declaration crashes (observed on macOS/Accelerate). When --enable-zdot or
# --disable-zdot is given explicitly we honor it; otherwise we autodetect the
# BLAS complex-return ABI and fall back to Fortran intrinsics (-D__NO_ZDOT)
# whenever the BLAS zdotc does not use the register-return convention.
# zdotu shares the same complex-return convention as zdotc, so probing zdotc
# alone is sufficient.
#
AC_ARG_ENABLE(zdot, [AS_HELP_STRING([--enable-zdot, --disable-zdot],
                           [Use (yes) or do not use (no) ZDOTC and ZDOTU from
                            BLAS. Default: autodetect the BLAS ABI.])],
                           [zdot_explicit=yes],
                           [zdot_explicit=no; enable_zdot=yes])
if test x"$zdot_explicit" = xyes; then
  if test x"$enable_zdot" != xyes; then
    FCFLAGS="$FCFLAGS -D__NO_ZDOT"
  fi
else
  AC_LANG_PUSH([Fortran])
  ac_save_LIBS="$LIBS"
  LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS"
  AC_MSG_CHECKING([whether BLAS zdotc uses the register-return ABI])
  AC_RUN_IFELSE([AC_LANG_SOURCE([[
      program conftest_zdotc
      implicit none
      complex(8) :: zx(2), zy(2), res
      complex(8) :: zdotc
      external zdotc
      zx(1) = (1d0, 0d0)
      zx(2) = (0d0, 1d0)
      zy(1) = (1d0, 0d0)
      zy(2) = (0d0, 1d0)
      res = zdotc(2, zx, 1, zy, 1)
      if (abs(dble(res) - 2d0) < 1d-8 .and. abs(aimag(res)) < 1d-8) then
        stop 0
      else
        stop 1
      end if
      end program conftest_zdotc
]])],
    [zdot_abi_ok=yes],
    [zdot_abi_ok=no],
    [zdot_abi_ok=cross])
  LIBS="$ac_save_LIBS"
  AC_LANG_POP([Fortran])
  if test x"$zdot_abi_ok" = xyes; then
    AC_MSG_RESULT([yes])
  elif test x"$zdot_abi_ok" = xcross; then
    AC_MSG_RESULT([assuming no (cross-compiling)])
    AC_MSG_WARN([cannot run a test program while cross-compiling; assuming the
                 BLAS zdotc/zdotu are unusable and defining -D__NO_ZDOT. Pass
                 --enable-zdot to force use of the BLAS routines.])
    FCFLAGS="$FCFLAGS -D__NO_ZDOT"
  else
    AC_MSG_RESULT([no])
    FCFLAGS="$FCFLAGS -D__NO_ZDOT"
  fi
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
