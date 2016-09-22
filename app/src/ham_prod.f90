MODULE ham_vals
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  nsite
  !
  REAL(8),SAVE :: &
  & Jx, &
  & Jy, &
  & Jz, &
  & Dz
  !
END MODULE ham_vals
!
!
!
MODULE ham_prod_mod
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC ham_prod
  !
CONTAINS
!
!
!
SUBROUTINE ham_prod(veci,veco)
  !
  USE shiftk_vals, ONLY : ndim, inham
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: veci(ndim)
  COMPLEX(8),INTENT(OUT) :: veco(ndim)
  !
  veco(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
  !
  IF(inham == "") THEN
     CALL ham_prod_onthefly(veci,veco)
  ELSE
     CALL ham_prod_compress(veci,veco)
  END IF
  !
END SUBROUTINE ham_prod
!
! Hamiltonian-vector product
!
SUBROUTINE ham_prod_compress(veci,veco)
  !
  USE shiftk_vals, ONLY : ndim, nham, ndiag, ham, ham_indx
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: veci(ndim)
  COMPLEX(8),INTENT(OUT) :: veco(ndim)
  !
  INTEGER :: iham
  !
  DO iham = 1, ndiag
     veco(ham_indx(1,iham)) = veco(ham_indx(1,iham)) &
     &          + ham(iham) * veci(ham_indx(2,iham))
  END DO
  !
  DO iham = ndiag + 1, nham
     veco(ham_indx(1,iham)) = veco(ham_indx(1,iham)) &
     &          + ham(iham) * veci(ham_indx(2,iham))
     veco(ham_indx(2,iham)) = veco(ham_indx(2,iham)) &
     &   + CONJG(ham(iham)) * veci(ham_indx(1,iham))
  END DO
  !
END SUBROUTINE ham_prod_compress
!
!
!
SUBROUTINE ham_prod_onthefly(veci,veco)
  !
  USE shiftk_vals, ONLY : ndim, threshold
  USE ham_vals, ONLY : Jx, Jy, Jz, Dz, nsite
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: veci(ndim)
  COMPLEX(8),INTENT(OUT) :: veco(ndim)
  !
  INTEGER :: isite, isite1, mask1, mask2, mask12, spin, idim
  COMPLEX(8) :: matrix
  !
  LOGICAL,ALLOCATABLE :: uu(:), ud(:), du(:), dd(:), para(:), anti(:)
  INTEGER,ALLOCATABLE :: pair(:)
  !
  ALLOCATE(uu(ndim), ud(ndim), du(ndim), dd(ndim), para(ndim), anti(ndim), pair(ndim))
  !
  DO isite = 1, nsite
     !
     isite1 = MOD(isite, nsite) + 1
!WRITE(*,*) isite, isite1
     !
     mask12 = 0
     mask12 = IBSET(mask12, isite  - 1)
     mask12 = IBSET(mask12, isite1 - 1)
     !
     mask1 = 0
     mask1 = IBSET(mask1, isite  - 1)
     !
     mask2 = 0
     mask2 = IBSET(mask2, isite1 - 1)
     !
     uu(1:ndim) = .FALSE.
     dd(1:ndim) = .FALSE.
     ud(1:ndim) = .FALSE.
     du(1:ndim) = .FALSE.
     para(1:ndim) = .FALSE.
     anti(1:ndim) = .FALSE.
     !
     DO idim = 1, ndim
        !
        spin = IAND(idim - 1, mask12)
        pair(idim) = IEOR(idim - 1, mask12) + 1
        !
        IF(spin == mask12) THEN
           uu(idim) = .TRUE.
           para(idim) = .TRUE.
        ELSE IF(spin == 0) THEN
           dd(idim) = .TRUE.
           para(idim) = .TRUE.
        ELSE IF(spin == mask1) THEN
           ud(idim) = .TRUE.
           anti(idim) = .TRUE.
        ELSE
           du(idim) = .TRUE.
           anti(idim) = .TRUE.
        END IF
        !
!WRITE(*,*) uu(idim), du(idim), ud(idim), dd(idim), para(idim), anti(idim), pair(idim)
     END DO
!stop
     !
     ! S_{i z} S_{i+1 z}
     !
     IF(ABS(Jz) > threshold) THEN
        !
        matrix = CMPLX(0.25d0 * Jz, 0d0, KIND(0d0))
        DO idim = 1, ndim
           IF(para(idim)) veco(idim) = veco(idim) + matrix * veci(idim)
        END DO
        !
        matrix = CMPLX(- 0.25d0 * Jz, 0d0, KIND(0d0))
        DO idim = 1, ndim
           IF(anti(idim)) veco(idim) = veco(idim) + matrix * veci(idim)
        END DO
        !
     END IF        
     !
     ! S_{i}^+ S_{i+1}^- + S_{i}^- S_{i+1}^+
     !
     IF(ABS(Jx + Jy) > threshold .OR. ABS(Dz) > threshold) THEN
        !
        matrix = CMPLX(0.25d0 * (Jx + Jy), 0.5d0 * Dz, KIND(0d0)) 
        DO idim = 1, ndim
           IF(du(idim)) veco(idim) = veco(idim) + matrix * veci(pair(idim))
        END DO
        !
        matrix = CMPLX(0.25d0 * (Jx + Jy), - 0.5d0 * Dz, KIND(0d0)) 
        DO idim = 1, ndim
           IF(ud(idim)) veco(idim) = veco(idim) + matrix * veci(pair(idim))
        END DO
        !
     END IF
     !
     ! S_{i}^+ S_{i+1}^+ + S_{i}^+ S_{i+1}^+
     !
     IF(ABS(Jx - Jy) > threshold) THEN
        !
        matrix = CMPLX(0.25d0 * (Jx - Jy), 0d0, KIND(0d0)) 
        DO idim = 1, ndim
           IF(para(idim)) veco(idim) = veco(idim) + matrix * veci(pair(idim))
        END DO
        !
     END IF
     !
  END DO
  !
  DEALLOCATE(uu, ud, du, dd, para, anti, pair)
  !
END SUBROUTINE ham_prod_onthefly
!
END MODULE ham_prod_mod
