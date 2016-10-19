!
!    Copyright 2016 Mitsuaki Kawamura
!
!    This file is part of ISSP Math Library.
!
!    ISSP Math Library is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ISSP Math Library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ISSP Math Library.  If not, see <http://www.gnu.org/licenses/>.
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
