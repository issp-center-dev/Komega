!
! ISSP Math Library - A library for solving linear systems in materials science
! Copyright (C) 2016 Mitsuaki Kawamura
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
! For more details, See `COPYING.LESSER' in the root directory of this library.
!
! Regression test for the komega_CG_R lz_conv lifecycle.
!
! lz_conv is ALLOCATEd unconditionally in komega_CG_R_init, so it must also be
! DEALLOCATEd unconditionally in komega_CG_R_finalize -- including the
! itermax == 0 case.  Before the fix, finalize only released lz_conv inside the
! IF(itermax > 0) block, so a second init with itermax == 0 aborted with
! "Attempting to allocate already allocated variable 'lz_conv'".
!
! This program performs two init/finalize cycles with itermax == 0.  Built with
! -fcheck=all it aborts on the unpatched library and runs to completion (exit 0)
! on the fixed one.
!
PROGRAM test_lzconv
  !
  USE komega_cg_r, ONLY : komega_CG_R_init, komega_CG_R_finalize
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: ndim = 4, nl = 4, nz = 2, itermax = 0
  REAL(8), PARAMETER :: threshold = 1d-8
  REAL(8) :: x(nl,nz), z(nz)
  INTEGER :: cyc
  !
  z(1) = 1d0
  z(2) = 2d0
  !
  DO cyc = 1, 2
     CALL komega_CG_R_init(ndim, nl, nz, x, z, itermax, threshold)
     CALL komega_CG_R_finalize()
  END DO
  !
  WRITE(*,'(a)') "test_lzconv: PASS (itermax==0 init/finalize x2 ok)"
  !
END PROGRAM test_lzconv
