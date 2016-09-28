MODULE ham_vals
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & ndiag,   & ! Diagonal components
  & nham,    & ! Non-zero elements of compressed Hamiltonian
  & nsite,   & ! Number of sites for the Heisenberg Chain
  & nsitep     ! Number of sites in inter-process region
  !
  REAL(8),SAVE :: &
  & Jx, & ! Heisenberg Jx
  & Jy, & ! Heisenberg Jy
  & Jz, & ! Heisenberg Jz
  & Dz    ! D-M interaction
  !
  LOGICAL,ALLOCATABLE,SAVE :: &
  & ud(:), & ! i=Up   & i+1=Down flag
  & du(:), & ! i=Down & i+1=Up   flag
  & para(:) ! i, i+1 Parallel flag
  !
  INTEGER,ALLOCATABLE,SAVE :: &
  & ham_indx(:,:), & ! row & column index of Hamiltonian
  & pair(:)          ! pair connected by S+S- or S+S+
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & ham(:) ! Compressed Hamiltonian
  !
END MODULE ham_vals
!
! Module contains routines for Hamiltonian-Vector product
!
MODULE ham_prod_mod
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC ham_prod, print_ham, onthefly_init, finalize_ham
  !
CONTAINS
!
! Driver routine for the Hamiltonian-Vector Product
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
! Hamiltonian-vector product with the Compressed representation
!
SUBROUTINE ham_prod_compress(veci,veco)
  !
  USE shiftk_vals, ONLY : ndim
  USE ham_vals, ONLY : nham, ndiag, ham, ham_indx
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
! Hamiltonian-vector product with On-The-Fly Hamiltonian generation
!
SUBROUTINE ham_prod_onthefly(veci,veco)
  !
  USE shiftk_vals, ONLY : ndim, almost0
  USE ham_vals, ONLY : Jx, Jy, Jz, Dz, nsite, ud, du, para, pair, nsitep
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: veci(ndim)
  COMPLEX(8),INTENT(OUT) :: veco(ndim)
  !
  INTEGER :: isite, isite1, mask1, mask2, mask12, spin, idim
  COMPLEX(8) :: matrix, cmatrix
  !
  ! Intra process region
  !
  DO isite = 1, nsite - nsitep
     !
     isite1 = MOD(isite, nsite) + 1
!WRITE(stdout,*) isite, isite1
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
     ud(1:ndim) = .FALSE.
     du(1:ndim) = .FALSE.
     para(1:ndim) = .FALSE.
     !
     DO idim = 1, ndim
        !
        spin = IAND(idim - 1, mask12)
        pair(idim) = IEOR(idim - 1, mask12) + 1
        !
        IF(spin == mask1) THEN
           ud(idim) = .TRUE.
        ELSE IF(spin == mask2) THEN
           du(idim) = .TRUE.
        ELSE
           para(idim) = .TRUE.
        END IF
        !
!WRITE(stdout,*) uu(idim), du(idim), ud(idim), dd(idim), para(idim), anti(idim), pair(idim)
     END DO
!stop
     !
     ! S_{i z} S_{i+1 z}
     !
     IF(ABS(Jz) > almost0) THEN
        !
        matrix = CMPLX(0.25d0 * Jz, 0d0, KIND(0d0))
        DO idim = 1, ndim
           IF(para(idim)) THEN
              veco(idim) = veco(idim) + matrix * veci(idim)
           ELSE
              veco(idim) = veco(idim) - matrix * veci(idim)
           END IF
        END DO
        !
     END IF        
     !
     ! S_{i}^+ S_{i+1}^- + S_{i}^- S_{i+1}^+
     !
     IF(ABS(Jx + Jy) > almost0 .OR. ABS(Dz) > almost0) THEN
        !
        matrix = CMPLX(0.25d0 * (Jx + Jy), 0.5d0 * Dz, KIND(0d0))
        cmatrix = CONJG(matrix)
        DO idim = 1, ndim
           IF(du(idim)) THEN
              veco(idim) = veco(idim) +  matrix * veci(pair(idim))
           ELSE IF(ud(idim)) THEN
              veco(idim) = veco(idim) + cmatrix * veci(pair(idim))
           END IF
        END DO
        !
     END IF
     !
     ! S_{i}^+ S_{i+1}^+ + S_{i}^- S_{i+1}^
     !
     IF(ABS(Jx - Jy) > almost0) THEN
        !
        matrix = CMPLX(0.25d0 * (Jx - Jy), 0d0, KIND(0d0)) 
        DO idim = 1, ndim
           IF(para(idim)) veco(idim) = veco(idim) + matrix * veci(pair(idim))
        END DO
        !
     END IF
     !
  END DO ! isite = 1, nsite
  !
  ! isite = nsite - nsitep .OR. isite = nsite
  !
#if defined(MPI)
  DO isite = 1, 2
     !
     IF(isite = 1) THEN
        mask1 = 0
        mask1 = IBSET(mask1, nsite - nsitep - 1)
        mask2 = 0
        mask2 = IBSET(mask2, 0)
     ELSE
        mask1 = 0
        mask1 = IBSET(mask1, 0)
        mask2 = 0
        mask2 = IBSET(mask2, nsitep - 1)
     END IF
     origin = IEOR(myrank, mask2)
     !
     CALL mpi_sendrecv(veci, ndim, MPI_DOUBLE_COMPLEX, isite1, 0, &
     &             veci_buf, ndim, MPI_DOUBLE_COMPLEX, isite1, 0, MPI_COMM_WORLD, status);
     !
     para(1:ndim) = .FALSE.
     !
     DO idim = 1, ndim
        !
        spin = IAND(idim - 1, mask1)
        pair(idim) = IEOR(idim - 1, mask1) + 1
        !
        IF(spin == mask1) para(idim) = .TRUE.
        !
     END DO
     !
     IF(IAND(myrank, mask2) == 0) THEN
        para(1:ndim) = .NOT. para(1:ndim)
     END IF
     !
     ! S_{i z} S_{i+1 z}
     !
     IF(ABS(Jz) > almost0) THEN
        !
        matrix = CMPLX(0.25d0 * Jz, 0d0, KIND(0d0))
        DO idim = 1, ndim
           IF(para(idim)) THEN
              veco(idim) = veco(idim) + matrix * veci(idim)
           ELSE
              veco(idim) = veco(idim) - matrix * veci(idim)
           END IF
        END DO
        !
     END IF        
     !
     ! S_{i}^+ S_{i+1}^- + S_{i}^- S_{i+1}^+
     !
     IF(ABS(Jx + Jy) > almost0 .OR. ABS(Dz) > almost0) THEN
        !
        IF(IAND(myrank, mask2) == mask2 .AND. isite == 1) THEN
           matrix = CMPLX(0.25d0 * (Jx + Jy), 0.5d0 * Dz, KIND(0d0))
        ELSE
           matrix = CMPLX(0.25d0 * (Jx + Jy), - 0.5d0 * Dz, KIND(0d0))
        END IF
        !
        DO idim = 1, ndim
           IF(.NOT. para(idim)) veco(idim) = veco(idim) + matrix * veci_buf(pair(idim))
        END DO
        !
     END IF
     !
     ! S_{i}^+ S_{i+1}^+ + S_{i}^- S_{i+1}^
     !
     IF(ABS(Jx - Jy) > almost0) THEN
        !
        matrix = CMPLX(0.25d0 * (Jx - Jy), 0d0, KIND(0d0)) 
        DO idim = 1, ndim
           IF(para(idim)) veco(idim) = veco(idim) + matrix * veci_buf(pair(idim))
        END DO
        !
     END IF
     !
  END DO ! isite = 1, 2
  !
  ! nsite - nsitep < isite < nsite
  !
  DO isite = 1, nsitep - 1
     !
     isite1 = isite + 1
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
     origin = IEOR(myrank, mask12)
     !
     CALL mpi_sendrecv(veci, ndim, MPI_DOUBLE_COMPLEX, isite1, 0, &
     &             veci_buf, ndim, MPI_DOUBLE_COMPLEX, isite1, 0, MPI_COMM_WORLD, status);
     !
     spin = IAND(myrank, mask12)
     !
     ! S_{i z} S_{i+1 z}
     !
     IF(ABS(Jz) > almost0) THEN
        !
        IF(spin == mask12 .OR. spin == 0) THEN
           matrix = CMPLX(0.25d0 * Jz, 0d0, KIND(0d0))
        ELSE
           matrix = CMPLX(- 0.25d0 * Jz, 0d0, KIND(0d0))
        END IF
        !
        DO idim = 1, ndim
           veco(idim) = veco(idim) + matrix * veci(idim)
        END DO
        !
     END IF        
     !
     ! S_{i}^+ S_{i+1}^- + S_{i}^- S_{i+1}^+
     !
     IF(ABS(Jx + Jy) > almost0 .OR. ABS(Dz) > almost0) THEN
        !
        IF(spin == mask1 .OR. spin == mask2) THEN
           !
           IF(spin == mask1) THEN
              matrix = CMPLX(0.25d0 * (Jx + Jy), - 0.5d0 * Dz, KIND(0d0))
           ELSE
              matrix = CMPLX(0.25d0 * (Jx + Jy),   0.5d0 * Dz, KIND(0d0))
           END IF
           !
           DO idim = 1, ndim
              veco(idim) = veco(idim) + matrix * veci_buf(idim)
           END DO
           !
        END IF
        !
     END IF
     !
     ! S_{i}^+ S_{i+1}^+ + S_{i}^- S_{i+1}^
     !
     IF(ABS(Jx - Jy) > almost0) THEN
        !
        IF(spin == mask12 .OR. spin == 0) THEN
           !
           matrix = CMPLX(0.25d0 * (Jx - Jy), 0d0, KIND(0d0)) 
           DO idim = 1, ndim
              IF(para(idim)) veco(idim) = veco(idim) + matrix * veci_buf(idim)
           END DO
           !
        END IF
        !
     END IF
     !
  END DO ! isite = 1, nsitep - 1
#endif
  !
END SUBROUTINE ham_prod_onthefly
!
! Allocate Flags and Pair index for On-The-Fly
!
SUBROUTINE onthefly_init()
  !
  USE shiftk_vals, ONLY : ndim
  USE ham_vals, ONLY : ud, du, para, pair
  !
  IMPLICIT NONE
  !
  ALLOCATE(ud(ndim), du(ndim), para(ndim), pair(ndim))
#if defined(MPI)
  ALLOCATE(veci_buf(ndim))
#endif
  !
END SUBROUTINE onthefly_init
!
! Deallocate Hamiltonian parameters
!
SUBROUTINE finalize_ham()
  !
  USE shiftk_vals, ONLY : inham
  USE ham_vals, ONLY : ud, du, para, pair, ham, ham_indx
  !
  IMPLICIT NONE
  !
  IF(inham == "") THEN
     DEALLOCATE(ud, du, para, pair)
#if defined(MPI)
     DEALLOCATE(veci_buf)
#endif
  ELSE
     DEALLOCATE(ham, ham_indx)
  END IF
  !
END SUBROUTINE finalize_ham
!
! Print Hamiltonian with the Matrix Market format
!
SUBROUTINE print_ham()
  !
  USE shiftk_vals, ONLY : ndim, almost0, myrank
  USE ham_vals, ONLY : nham
  !
  IMPLICIT NONE
  !
  INTEGER :: idim, jdim, fo = 21
  COMPLEX(8),ALLOCATABLE :: veci(:), veco(:)
  !
  ALLOCATE(veci(ndim), veco(ndim))
  !
  IF(myrank == 0) THEN
     OPEN(fo, file = "zvo_Ham.dat")
     WRITE(fo,*) "%%MatrixMarket matrix coordinate complex hermitian"
  END IF
  !
  nham = 0
  DO idim = 1, ndim
     !
     veci(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
     veci(  idim) = CMPLX(1d0, 0d0, KIND(0d0))
     !
     CALL ham_prod(veci, veco)
     !
     nham = nham + COUNT(ABS(veco(idim:ndim)) > almost0)
     !
  END DO
  !
  IF(myrank == 0) WRITE(fo,'(3i10)') ndim, ndim, nham 
  !
  DO idim = 1, ndim
     !
     veci(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
     veci(  idim) = CMPLX(1d0, 0d0, KIND(0d0))
     !
     CALL ham_prod(veci, veco)
     !
     IF(myrank == 0) THEN
        DO jdim = idim, ndim
           IF(ABS(veco(jdim)) > almost0) &
           &  WRITE(fo,'(2i10,2f15.8)') jdim, idim, DBLE(veco(jdim)), AIMAG(veco(jdim))
        END DO
     END IF
     !
  END DO
  !
  IF(myrank == 0) THEN
     CLOSE(fo)
  END IF
  !
  DEALLOCATE(veci, veco)
  !
END SUBROUTINE print_ham
!
END MODULE ham_prod_mod
