#!/bin/sh
#
# ISSP Math Library - A library for solving linear systems in materials science
# Copyright (C) 2016 Mitsuaki Kawamura
#
# Regression test runner, invoked by "make check" (check-local in test/Makefile).
#
# It performs two kinds of checks:
#
#   1. Numerical correctness of the four shifted-Krylov drivers
#      (solve_rr/cr/rc/cc.x).  Each is run on a small, well-conditioned system
#      with a generous iteration budget (regression_*.in) and must
#        (a) exit 0,
#        (b) report "Converged in iteration",
#        (c) produce a final residual ||(z-H) x - b|| below RESIDUAL_TOL.
#      The drivers build their matrix with RANDOM_NUMBER (no fixed seed, so the
#      Green's function values are not reproducible across platforms or runs);
#      the residual ||(z-H) x - b||, however, must be ~0 for ANY well-conditioned
#      system once the solver converges, so it is the portable, deterministic
#      correctness signal -- that is what we assert against a tolerance, rather
#      than a platform-dependent stored Green's function.
#      With itermax=50 / threshold=1d-12 on these small (ndim=5) systems the
#      solver converges to machine precision: the worst residual observed over
#      hundreds of runs is ~1e-12, so the 1e-8 tolerance keeps a >=4 orders of
#      magnitude safety margin against the random matrix.
#      The random Hamiltonian is positive semi-definite (A^T A / A^H A, real
#      eigenvalues >= 0).  For the real-frequency drivers the shifts in
#      regression_real.in are therefore all negative, so (z - H) can never be
#      singular (no shift collides with the spectrum); the complex-frequency
#      inputs already keep Im(z) = 1, which bounds the conditioning.  This is
#      what makes convergence reliable regardless of the (unseeded) matrix.
#
#   2. The komega_CG_R lz_conv itermax==0 lifecycle regression (test_lzconv),
#      compiled with -fcheck=all and required to run to completion.
#
# Exit status is non-zero if any check fails, so "make check" fails the build.
#
# Variables are normally provided by the Makefile, with fall-backs for a manual
# run from the test build directory.
#
set -u

srcdir=${srcdir:-.}
top_builddir=${top_builddir:-..}
FC=${FC:-gfortran}
FCFLAGS=${FCFLAGS:-}
LAPACK_LIBS=${LAPACK_LIBS:-}
BLAS_LIBS=${BLAS_LIBS:-}
LIBTOOL=${LIBTOOL:-${top_builddir}/libtool}
RESIDUAL_TOL=${RESIDUAL_TOL:-1e-8}

fail=0

note() { printf '%s\n' "$*"; }

# Largest absolute value in the "Residual Vector" block of a driver's stdout.
max_residual() {
  awk '
    /Residual Vector/ {f=1; next}
    /#####/           {f=0}
    f { for (i=1;i<=NF;i++){ v=$i+0; if (v<0) v=-v; if (v>m) m=v } }
    END { printf "%.6e", m+0 }
  '
}

run_driver() {
  # $1 = driver basename, $2 = input file (in srcdir)
  drv=$1
  inp="${srcdir}/$2"
  out=$(./"${drv}".x < "${inp}" 2>&1)
  rc=$?
  if [ $rc -ne 0 ]; then
    note "FAIL ${drv}: nonzero exit (${rc})"
    note "${out}" | tail -n 5
    fail=1
    return
  fi
  if ! printf '%s\n' "${out}" | grep -q "Converged in iteration"; then
    note "FAIL ${drv}: did not converge"
    note "${out}" | tail -n 5
    fail=1
    return
  fi
  res=$(printf '%s\n' "${out}" | max_residual)
  ok=$(awk -v r="${res}" -v t="${RESIDUAL_TOL}" 'BEGIN{ print (r+0 <= t+0) ? "yes" : "no" }')
  if [ "${ok}" != "yes" ]; then
    note "FAIL ${drv}: residual ${res} > tol ${RESIDUAL_TOL}"
    fail=1
    return
  fi
  note "PASS ${drv}: converged, max residual ${res} <= ${RESIDUAL_TOL}"
}

note "===== Komega regression tests ====="

# 1. Numerical correctness of the four drivers.
rm -f restart.dat
run_driver solve_rr regression_real.in
rm -f restart.dat
run_driver solve_cr regression_real.in
rm -f restart.dat
run_driver solve_rc regression_complex.in
rm -f restart.dat
run_driver solve_cc regression_complex.in
rm -f restart.dat

# 2. lz_conv itermax==0 lifecycle regression (built with -fcheck=all).
note "----- komega_CG_R lz_conv itermax==0 regression -----"
if "${LIBTOOL}" --mode=link --tag=FC ${FC} ${FCFLAGS} -fcheck=all -fbacktrace \
     -I"${top_builddir}/src" -o test_lzconv.x "${srcdir}/test_lzconv.F90" \
     "${top_builddir}/src/libkomega.la" ${LAPACK_LIBS} ${BLAS_LIBS} \
     >test_lzconv_build.log 2>&1; then
  if "${LIBTOOL}" --mode=execute ./test_lzconv.x; then
    note "PASS test_lzconv"
  else
    note "FAIL test_lzconv: aborted (lz_conv not deallocated for itermax==0?)"
    fail=1
  fi
else
  note "FAIL test_lzconv: build failed"
  cat test_lzconv_build.log
  fail=1
fi
rm -f test_lzconv.x test_lzconv_build.log

note "==================================="
if [ $fail -ne 0 ]; then
  note "RESULT: FAIL"
  exit 1
fi
note "RESULT: PASS"
exit 0
