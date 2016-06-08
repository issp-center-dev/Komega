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
