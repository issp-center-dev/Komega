#!/bin/sh
#
# ISSP Math Library - A library for solving linear systems in materials science
# Copyright (C) 2016 Mitsuaki Kawamura
#
# Smoke test for the configure-time zdot (BLAS complex-return ABI) handling.
#
# It exports a pristine copy of the repository (git archive) into a temporary
# directory and runs ./configure three times, asserting:
#
#   * --disable-zdot  : autodetection is skipped, -D__NO_ZDOT IS defined, the
#                       build succeeds, and "make check" passes (the Fortran
#                       intrinsic fallback works with any BLAS, so this gives
#                       runtime coverage of the -D__NO_ZDOT path).
#   * --enable-zdot   : autodetection is skipped, -D__NO_ZDOT is NOT defined,
#                       and the build succeeds.  The complex-dot drivers are
#                       only valid where the BLAS uses the register-return ABI,
#                       so "make check" is additionally run for this mode only
#                       when KOMEGA_RUN_ENABLE=1 (set by the Linux CI job);
#                       otherwise this mode asserts the configure decision +
#                       build only.
#   * bare configure  : autodetection runs (the "checking whether BLAS zdotc
#                       uses the register-return ABI" line appears), the build
#                       succeeds, and -- because the bare default always picks
#                       the safe/correct path -- "make check" passes.
#
# The expected bare-mode autodetect outcome is platform dependent (a
# register-return BLAS such as OpenBLAS/MKL -> no -D__NO_ZDOT; the f2c-ABI
# standard BLAS / Accelerate on macOS -> -D__NO_ZDOT).  Set the environment
# variable KOMEGA_EXPECT_BARE to "no" or "yes" to assert a specific outcome
# (recommended on CI, where the BLAS ABI is known); the default "auto" accepts
# either, so the build still has to succeed but the decision is not pinned.
#
# This script must be run from inside the git work tree (CI checks out the repo
# before invoking it).  It does not touch the live build tree.
#
set -u

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "check_configure_zdot: must be run inside the Komega git work tree" >&2
  exit 2
fi
srcroot=${srcroot:-$(git rev-parse --show-toplevel)}

FC=${FC:-gfortran}
MAKE=${MAKE:-make}
KOMEGA_EXPECT_BARE=${KOMEGA_EXPECT_BARE:-auto}
case "${KOMEGA_EXPECT_BARE}" in
  yes|no|auto) ;;
  *) echo "check_configure_zdot: KOMEGA_EXPECT_BARE must be yes, no or auto (got '${KOMEGA_EXPECT_BARE}')" >&2
     exit 2 ;;
esac
work=$(mktemp -d 2>/dev/null || mktemp -d -t komega_zdot)
trap 'rm -rf "${work}"' 0

fail=0

# Export a pristine tree for one configure run.
export_tree() {
  dest=$1
  mkdir -p "${dest}"
  ( cd "${srcroot}" && git archive HEAD ) | tar -x -C "${dest}"
}

# Does the configured tree define -D__NO_ZDOT?  Read the generated src/Makefile.
# Use "grep -e" (POSIX) so the leading '-' in the pattern is not parsed as an
# option.
has_no_zdot() {
  grep -q -e '-D__NO_ZDOT' "$1/src/Makefile"
}

run_mode() {
  # $1 = label, $2 = configure flag ("" for bare), $3 = expect macro (yes/no/auto)
  label=$1; flag=$2; expect=$3
  dir="${work}/${label}"
  export_tree "${dir}"
  echo "----- mode: ${label} (configure ${flag:-<bare>}) -----"

  if ! ( cd "${dir}" && ./configure FC="${FC}" ${flag} ) >"${dir}/configure.log" 2>&1; then
    echo "FAIL ${label}: configure failed"; tail -n 15 "${dir}/configure.log"; fail=1; return
  fi

  if [ "${label}" = "bare" ]; then
    if ! grep -q "whether BLAS zdotc uses the register-return ABI" "${dir}/configure.log"; then
      echo "FAIL ${label}: autodetection check did not run"; fail=1; return
    fi
    echo "  autodetect: $(grep 'whether BLAS zdotc' "${dir}/configure.log" | sed 's/^.*\.\.\. //')"
  fi

  if has_no_zdot "${dir}"; then macro=yes; else macro=no; fi
  case "${expect}" in
    yes) if [ "${macro}" != yes ]; then echo "FAIL ${label}: expected -D__NO_ZDOT, absent"; fail=1; return; fi ;;
    no)  if [ "${macro}" != no  ]; then echo "FAIL ${label}: -D__NO_ZDOT unexpectedly defined"; fail=1; return; fi ;;
    auto) : ;;  # bare: decision is platform-dependent, accept either
  esac
  echo "  -D__NO_ZDOT defined: ${macro} (expected: ${expect})"

  if ! ( cd "${dir}" && ${MAKE} ) >"${dir}/make.log" 2>&1; then
    echo "FAIL ${label}: build failed"; tail -n 15 "${dir}/make.log"; fail=1; return
  fi
  echo "  build: ok"

  # Run the drivers end-to-end where the BLAS ABI permits.  The bare default
  # always picks a working configuration and the --disable-zdot fallback uses
  # the Fortran intrinsics (DOT_PRODUCT/SUM), which work with any BLAS, so
  # "make check" is run for both -- this gives runtime coverage of the
  # -D__NO_ZDOT fallback path.  The forced --enable-zdot path only works when
  # the BLAS uses the register-return ABI (e.g. OpenBLAS/MKL on Linux), so it is
  # exercised only when the caller opts in via KOMEGA_RUN_ENABLE=1 (Linux CI).
  run_check=no
  [ "${label}" = "bare" ] && run_check=yes
  [ "${label}" = "disable" ] && run_check=yes
  [ "${label}" = "enable" ] && [ "${KOMEGA_RUN_ENABLE:-0}" = "1" ] && run_check=yes
  if [ "${run_check}" = yes ]; then
    if ! ( cd "${dir}" && ${MAKE} check ) >"${dir}/check.log" 2>&1; then
      echo "FAIL ${label}: make check failed"; tail -n 20 "${dir}/check.log"; fail=1; return
    fi
    echo "  make check: pass"
  fi

  echo "PASS ${label}"
}

echo "===== configure zdot-autodetect smoke test ====="
echo "(bare-mode expected autodetect outcome: ${KOMEGA_EXPECT_BARE})"
run_mode disable "--disable-zdot" yes
run_mode enable  "--enable-zdot"  no
run_mode bare    ""               "${KOMEGA_EXPECT_BARE}"
echo "==============================================="

if [ ${fail} -ne 0 ]; then
  echo "RESULT: FAIL"; exit 1
fi
echo "RESULT: PASS"; exit 0
