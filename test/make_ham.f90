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
MODULE make_ham
  !
  IMPLICIT NONE
  !
CONTAINS
!
SUBROUTINE make_random_pd(a,n)
  !
  USE mathlib, ONLY : dgemm, dcopy
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n
  REAL(8),INTENT(INOUT) :: a(n,n)
  !
  REAL(8),ALLOCATABLE :: a0(:,:)
  !
  ALLOCATE(a0(n,n))
  !
  CALL RANDOM_SEED()
  !
  CALL RANDOM_NUMBER(a(1:n,1:n))
  !
  CALL dgemm("T", "N", n, n, n, 1d0, a, n, a, n, 0d0, a0, n)
  !
  CALL dcopy(n*n,a0,1,a,1)
  !
  DEALLOCATE(a0)
  !
END SUBROUTINE make_random_pd
!
END MODULE make_ham
