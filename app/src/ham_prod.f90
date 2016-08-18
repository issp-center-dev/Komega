MODULE ham_prod
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Hamiltonian-vector product
!
SUBROUTINE ham_prod_compress(veci,veco)
  !
  USE shiftk_vals, ONLY : ndim, nham, ham, ham_indx
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: veci(ndim)
  COMPLEX(8),INTENT(OUT) :: veco(ndim)
  !
  INTEGER :: iham
  !
  veco(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
  !
  DO iham = 1, nham
     veco(ham_indx(1,iham)) = veco(ham_indx(1,iham)) &
     &          + ham(iham) * veci(ham_indx(2,iham))
  END DO
  !
END SUBROUTINE ham_prod_compress
!
END MODULE ham_prod
