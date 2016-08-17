MODULE shiftk_vals
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & nham,    & ! Non-zero elements of compressed Hamiltonian
  & ndim,    & ! Size of Hilvert space
  & nl,      & ! Dimention of x
  & nomega,      & ! Number of frequencies
  & maxloops, & ! Max. number of iteraction
  & iter_old   ! Number of iteraction of previous run
  !
  REAL(8),SAVE :: &
  & threshold ! Convergence Threshold
  !
  COMPLEX(8),SAVE :: &
  & z_seed ! Seed frequency
  !
  INTEGER,ALLOCATABLE,SAVE :: &
  & ham_indx(:,:) ! row & column index of Hamiltonian
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & z(:)         ! (nomega): Frequency
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & ham(:), & ! Compressed Hamiltonian
  & rhs(:), &
  & v12(:), v2(:), & ! (ndim): Working vector
  & v14(:), v4(:), & ! (ndim): Working vector
  & r_l(:), & ! (nl) : Projeccted residual vector 
  & x(:,:) ! (nl,nomega) : Projected result 
  !
  ! Variables for Restart
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & alpha(:), beta(:) ! (iter_old) 
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & r_l_save(:,:) ! (nl,iter_old) Projected residual vectors
  !
  CHARACTER(256),SAVE :: &
  & inham, &
  & invec, &
  & calctype
  !
  LOGICAL,SAVE :: &
  lBiCG, &
  outrestart
  !
END MODULE shiftk_vals
!
! Routines
!
MODULE shiftk_routines
  !
  IMPLICIT NONE
  !
CONTAINS
  !
SUBROUTINE input_filename()
  !
  USE shiftk_vals, ONLY : inham, invec
  !
  IMPLICIT NONE
  !
  NAMELIST /filename/ inham, invec
  !
  inham = "zvo_Ham.dat"
  invec = "zvo_Excited.dat"
  !
  READ(*,filename,err=100)
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input FileName  ##########"
  WRITE(*,*)
  WRITE(*,*) "  Hamiltonian file : ", TRIM(inham)
  WRITE(*,*) "  Excited state file : ", TRIM(invec)
  !
  return
  !
100 write(*,*) "Stop in INPUT_FILENAME. reading namelist FILENAME"
  !
  stop
  !
END SUBROUTINE input_filename
!
!
!
SUBROUTINE input_parameter()
  !
  USE shiftk_vals, ONLY : nomega, maxloops, outrestart, threshold, ndim, calctype, z
  !
  IMPLICIT NONE
  !
  INTEGER :: convfactor, iomega
  COMPLEX(8) :: omegamax, omegamin
  NAMELIST /parameter/ omegamax, omegamin, nomega, maxloops, calctype, convfactor, outrestart
  !
  maxloops = ndim
  calctype = "normal"
  convfactor = 8
  nomega = 10
  omegamin = CMPLX(0d0, 1d0, KIND(0d0))
  omegamax = CMPLX(1d0, 1d0, KIND(0d0))
  outrestart = .FALSE.
  !
  READ(*,parameter,err=100)
  threshold = 10d0**(-convfactor)
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Parameter  ##########"
  WRITE(*,*)
  WRITE(*,*) "           Max. of Omega : ", omegamax
  WRITE(*,*) "           Min. of Omega : ", omegamin
  WRITE(*,*) "           Num. of Omega : ", nomega
  WRITE(*,*) "  Maximum number of loop : ", maxloops
  WRITE(*,*) "   Convergence Threshold : ", threshold
  WRITE(*,'(a,a)') "         Calculation type : ", calctype
  !
  ALLOCATE(z(nomega))
  z(1) = omegamin
  DO iomega = 2, nomega
     z(iomega) = omegamin + (omegamax - omegamin) * DBLE(iomega - 1) / DBLE(nomega - 1)
  END DO
  !
  return
  !
100 write(*,*) "Stop in INPUT_PARAMETER. reading namelist PARAMETER"
  !
  stop
  !
END SUBROUTINE input_parameter
!
! Input Hamiltonian with the Matrix Market Format
!
SUBROUTINE input_hamiltonian()
  !
  USE shiftk_vals, ONLY : ndim, ham, ham_indx, nham, inham, lBiCG
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ndim2, iham
  REAL(8) :: ham_r, ham_i
  CHARACTER(100) :: ctmp
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Hamiltonian  ##########"
  WRITE(*,*)
  !
  OPEN(fi, file = TRIM(inham))
  READ(fi, *) ctmp
  READ(fi,*) ndim, ndim2, nham
  WRITE(*,*) "          Dim. of Hamiltonian : ", ndim, ndim2
  WRITE(*,*) "  Num. of Non-Zero Components : ", nham
  !
  IF(ndim2 /= ndim) THEN
     WRITE(*,*) "ERROR ! Hamiltonian is not square."
     STOP
  END IF
  !
  ALLOCATE(ham(nham), ham_indx(2,nham))
  !
  DO iham = 1, nham
     READ(fi,*) ham_indx(1,iham), ham_indx(2,iham), ham_r, ham_i
     ham(iham) = CMPLX(ham_r, ham_i, KIND(0d0))
  END DO
  !
  CLOSE(fi)
  !
  IF(MAXVAL(ABS(AIMAG(ham(1:nham)))) > 1d-8) THEN
     WRITE(*,*) "  BiCG mathod is used."
     lBiCG = .TRUE.
  ELSE
     WRITE(*,*) "  COCG mathod is used."
     lBiCG = .FALSE.
  END IF
  !
END SUBROUTINE input_hamiltonian
!
! Input Excited State
!
SUBROUTINE input_rhs_vector()
  !
  USE shiftk_vals, ONLY : ndim, rhs, invec
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ndim2, idim
  REAL(8) :: rhs_r, rhs_i
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Right Hand Side Vector ##########"
  WRITE(*,*)
  OPEN(fi, file = TRIM(invec))
  READ(fi,*) ndim2
  WRITE(*,*) "  Dim. of RHS vector : ", ndim2
  !
  IF(ndim2 /= ndim) THEN
     WRITE(*,*) "ERROR ! Dimension is Incorrect."
     STOP
  END IF
  !
  ALLOCATE(rhs(ndim))
  !
  DO idim = 1, ndim
     READ(fi,*) rhs_r, rhs_i
     rhs(idim) = CMPLX(rhs_r, rhs_i, KIND(0d0))
  END DO
  !
  CLOSE(fi)
  !
END SUBROUTINE input_rhs_vector
!
! Input restart variables from file
!
SUBROUTINE input_restart_parameter()
  !
  USE shiftk_vals, ONLY : iter_old, alpha, beta, z_seed, nl, r_l_save
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, iter, il
  REAL(8) :: z_seed_r, z_seed_i, alpha_r, alpha_i, beta_r, beta_i
  REAL(8) :: r_l_save_r, r_l_save_i
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Restart Parameter  ##########"
  WRITE(*,*)
  !
  OPEN(fi, file = 'TriDiagComp.dat')
  !
  READ(fi,*) iter_old
  WRITE(*,*) "  Num. of Iteration (Previous Run) : ", iter_old
  !
  ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl, iter_old))
  !
  READ(fi,*) z_seed_r, z_seed_i
  WRITE(*,*) "  Previous Omega_Seed : ", z_seed_r, z_seed_i
  z_seed = CMPLX(z_seed_r, z_seed_i, KIND(0d0))
  !
  DO iter = 1, iter_old
     !
     READ(fi,*) alpha_r, alpha_i, beta_r, beta_i
     alpha(iter) = CMPLX(alpha_r, alpha_i, KIND(0d0))
     beta(iter) = CMPLX(beta_r, beta_i, KIND(0d0))
     !
  END DO
  !
  DO iter = 1, iter_old
     !
     DO il = 1, nl
        READ(fi,*) r_l_save_r, r_l_save_i
        r_l_save(il,iter) = CMPLX(r_l_save_r, r_l_save_i, KIND(0d0))
     END DO
     !
  END DO
  !
  CLOSE(fi)
  !
END SUBROUTINE input_restart_parameter
!
! Input restart variables from file
!
SUBROUTINE input_restart_vector()
  !
  USE shiftk_vals, ONLY : ndim, v2, v12, v4, v14, lBiCG
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ndim2, idim
  REAL(8) :: v2_r, v2_i, v12_r, v12_i
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Restart Vector  ##########"
  WRITE(*,*)
  !
  OPEN(fi, file = 'ResVec.dat')
  !
  READ(fi,*) ndim2
  WRITE(*,*) "  Dim. of Residual vector : ", ndim2
  !
  IF(ndim2 /= ndim) THEN
     WRITE(*,*) "ERROR ! Dimension is Incorrect."
     STOP
  END IF
  !
  DO idim = 1, ndim
     READ(fi,*) v2_r, v2_i, v12_r, v12_i
     v2(idim) = CMPLX(v2_r, v2_i, KIND(0d0))
     v12(idim) = CMPLX(v12_r, v12_i, KIND(0d0))
  END DO
  !
  IF(lBiCG) THEN
     READ(fi,*) v2_r, v2_i, v12_r, v12_i
     v4(idim) = CMPLX(v2_r, v2_i, KIND(0d0))
     v14(idim) = CMPLX(v12_r, v12_i, KIND(0d0))
  END IF
  !
  CLOSE(fi)
  !
END SUBROUTINE input_restart_vector
!
! Hamiltonian-vector product
!
SUBROUTINE ham_prod(veci,veco)
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
END SUBROUTINE ham_prod
!
! Input restart variables from file
!
SUBROUTINE output_restart_parameter()
  !
  USE shiftk_vals, ONLY : iter_old, alpha, beta, z_seed, r_l_save, nl
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, iter, il
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Output Restart Parameter  ##########"
  WRITE(*,*)
  !
  OPEN(fo, file = 'TriDiagComp.dat')
  !
  WRITE(fo,*) iter_old
  WRITE(*,*) "  Num. of Iteration (Current Run) : ", iter_old
  !
  WRITE(*,'(a,2e13.5)') "   Current Omega_Seed : ", z_seed
  WRITE(fo,'(2e25.16)') z_seed
  !
  DO iter = 1, iter_old
     !
     WRITE(fo,'(4e25.16)') alpha(iter), beta(iter)
     !
  END DO
  !
  DO iter = 1, iter_old
     !
     DO il = 1, nl
        WRITE(fo,'(2e25.16)') r_l_save(il,iter)
     END DO
     !
  END DO
  !
  CLOSE(fo)
  !
END SUBROUTINE output_restart_parameter
!
! Input restart variables from file
!
SUBROUTINE output_restart_vector()
  !
  USE shiftk_vals, ONLY : ndim, v2, v12, v4, v14, lBiCG
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, idim
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Output Restart Vector  ##########"
  WRITE(*,*)
  !
  OPEN(fo, file = 'ResVec.dat')
  !
  WRITE(fo,*) ndim
  WRITE(*,*) "  Dim. of Residual vector : ", ndim
  !
  DO idim = 1, ndim
     WRITE(fo,'(4e25.16)') v2(idim), v12(idim)
  END DO
  !
  IF(lBiCG) THEN
     DO idim = 1, ndim
        WRITE(fo,'(4e25.16)') v4(idim), v14(idim)
     END DO
  END IF
  !
  CLOSE(fo)
  !
END SUBROUTINE output_restart_vector
!
! Check Result
!
SUBROUTINE output_result()
  !
  USE shiftk_vals, ONLY : x, z, nomega
  !
  IMPLICIT NONE
  !
  INTEGER :: iz, fo = 20
  !
  OPEN(fo, file = "dynamicalG.dat")
  !
  DO iz = 1, nomega
     !
     write(fo,'(4e13.5)') DBLE(z(iz)), AIMAG(z(iz)), DBLE(x(1,iz)), AIMAG(x(1,iz))
     !
  END DO
  !
  CLOSE(fo)
  !
END SUBROUTINE output_result
!
! Check Result
!
SUBROUTINE output_result_debug()
  !
  USE shiftk_vals, ONLY : v2, ndim, x, rhs, z, nomega
  !
  IMPLICIT NONE
  !
  INTEGER :: iz, fo = 20
  CHARACTER(100) :: cndim, form
  COMPLEX(8) :: Gii
  !
  WRITE(cndim,*) ndim * 2
  WRITE(form,'(a,a,a)') "(", TRIM(ADJUSTL(cndim)), "e13.5)"
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Check Results  #####"
  WRITE(*,*)
  !
  WRITE(*,*) "  Residual Vector"
  !
  DO iz = 1, nomega
     !
     CALL ham_prod(x(1:ndim,iz), v2)
     v2(1:ndim) = z(iz) * x(1:ndim,iz) - v2(1:ndim) - rhs(1:ndim)
     !
     !write(*,form) v2(1:ndim)
     write(*,'(a,i5,a,2e13.5,a,e13.5)') "DEBUG (", iz, "), omega = ", z(iz), &
     &            ", Res. = ", dble(dot_product(v2,v2))
     !
  END DO
  !
  OPEN(fo, file = "dynamicalG.dat")
  !
  DO iz = 1, nomega
     !
     Gii = dot_product(rhs,x(1:ndim,iz))
     write(fo,'(4e13.5)') DBLE(z(iz)), AIMAG(z(iz)), DBLE(Gii), AIMAG(Gii)
     !
  END DO
  !
  CLOSE(fo)
  !
END SUBROUTINE output_result_debug
!
END MODULE shiftk_routines
!
!
!
PROGRAM shiftk
  !
  USE shifted_cocg, ONLY : COCG_init, COCG_restart, COCG_update, &
  &                        COCG_getcoef, COCG_getvec, COCG_finalize
  USE shifted_bicg, ONLY : BiCG_init, BiCG_restart, BiCG_update, &
  &                        BiCG_getcoef, BiCG_getvec, BiCG_finalize
  USE shiftk_routines, ONLY : input_filename, input_hamiltonian, input_rhs_vector, ham_prod, &
  &                           input_parameter, input_restart_parameter, input_restart_vector, &
  &                           output_result, output_result_debug, output_restart_parameter, output_restart_vector
  USE shiftk_vals, ONLY : alpha, beta, calctype, ndim, nomega, maxloops, iter_old, lBiCG, nl, &
  &                       outrestart, rhs, v12, v2, v14, v4, rhs, r_l, r_l_save, threshold, x, z, z_seed
  !
  IMPLICIT NONE
  !
  ! Variables for Restart
  !
  INTEGER :: &
  & iter,    & ! Counter for Iteration
  & status(3)
  !
  ! Input Filename
  !
  CALL input_filename()
  CALL input_hamiltonian()
  CALL input_rhs_vector()
  !
  nl = 1
  !nl = ndim
  CALL input_parameter()
  !
  ALLOCATE(v12(ndim), v2(ndim), r_l(nl), x(nl,nomega))
  IF(lBiCG) ALLOCATE(v14(ndim), v4(ndim))
  !
  IF(TRIM(calctype) == "recalc" .OR. TRIM(calctype) == "restart") THEN
     !
     CALL input_restart_parameter()
     IF(TRIM(calctype) == "restart") CALL input_restart_vector()
     maxloops = MAX(maxloops, iter_old)
     !
     WRITE(*,*)
     WRITE(*,*) "##########  CG Restart  ##########"
     WRITE(*,*)
     !
     IF(outrestart == .TRUE.) THEN
        IF(lBiCG) THEN
           CALL BiCG_restart(ndim, nl, nomega, x, z, maxloops, threshold, status, &
           &                 iter_old, v2, v12, v4, v14, alpha, beta, z_seed, r_l_save)
        ELSE
           CALL COCG_restart(ndim, nl, nomega, x, z, maxloops, threshold, status, &
           &                 iter_old, v2, v12,          alpha, beta, z_seed, r_l_save)
        END IF
     ELSE
        IF(lBiCG) THEN
           CALL BiCG_restart(ndim, nl, nomega, x, z, 0,        threshold, status, &
           &                 iter_old, v2, v12, v4, v14, alpha, beta, z_seed, r_l_save)
        ELSE
           CALL COCG_restart(ndim, nl, nomega, x, z, 0,        threshold, status, &
           &                 iter_old, v2, v12,          alpha, beta, z_seed, r_l_save)
        END IF
     END IF
     DEALLOCATE(alpha, beta, r_l_save)
     !
     IF(iter_old == maxloops .OR. TRIM(calctype) == "recalc") GOTO 10
     !
  ELSE IF(TRIM(calctype) == "normal") THEN
     !
     WRITE(*,*)
     WRITE(*,*) "##########  CG Initialization  ##########"
     WRITE(*,*)
     !
     v2(1:ndim) = rhs(1:ndim)
     !
     IF(outrestart == .TRUE.) THEN
        IF(lBiCG) THEN
           CALL BiCG_init(ndim, nl, nomega, x, z, maxloops, threshold, status)
        ELSE
           CALL COCG_init(ndim, nl, nomega, x, z, maxloops, threshold, status)
        END IF
     ELSE
        IF(lBiCG) THEN
           CALL BiCG_init(ndim, nl, nomega, x, z, 0,        threshold, status)
        ELSE
           CALL COCG_init(ndim, nl, nomega, x, z, 0,        threshold, status)
        END IF
     END IF
     !
  ELSE
     !
     WRITE(*,*) "ERROR ! calctype = ", TRIM(calctype)
     STOP
     !
  END IF
  !
  ! COCG/BiCG Loop
  !
  WRITE(*,*)
  WRITE(*,*) "#####  BiCG Iteration  #####"
  WRITE(*,*)
  !
  DO iter = 1, maxloops
     !
     ! Projection of Residual vector into the space
     ! spaned by left vectors
     !
     IF(ndim == nl) THEN
        r_l(1:ndim) = v2(1:ndim)
     ELSE
        r_l(1) = DOT_PRODUCT(rhs, v2)
     END IF
     !
     ! Matrix-vector product
     !
     CALL ham_prod(v2, v12)
     IF(lBiCG) CALL ham_prod(v4, v14)
     !
     ! Update result x with COCG
     !
     IF(lBiCG) THEN
        CALL BiCG_update(v12, v2, v14, v4, x, r_l, status)
     ELSE
        CALL COCG_update(v12, v2,          x, r_l, status)
     END IF
     !
     WRITE(*,'(a,4i,e13.5)') "  DEBUG : ", iter, status, DBLE(v12(1))
     IF(status(1) /= 0) EXIT
     !
  END DO
  !
  IF(status(1) > 0) THEN
     WRITE(*,*) "  Converged in iteration ", status(1)
  ELSE
     WRITE(*,*) "  Not Converged in iteration ", -status(1)
  END IF
  iter_old = abs(status(1))
  !
10 CONTINUE
  !
  ! Get these vectors for restart in the Next run
  !
  IF(outrestart == .TRUE.) THEN
     !
     ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl, iter_old))
     !
     IF(lBiCG) THEN
        CALL BiCG_getcoef(alpha, beta, z_seed, r_l_save)
        CALL BiCG_getvec(v12,v14)
     ELSE
        CALL COCG_getcoef(alpha, beta, z_seed, r_l_save)
        CALL COCG_getvec(v12)
     END IF
     !
     CALL output_restart_parameter()
     CALL output_restart_vector()
     !
     DEALLOCATE(alpha, beta, r_l_save)
     !     
  END IF
  !
  ! Deallocate all intrinsic vectors
  !
  IF(lBiCG) THEN
     CALL BiCG_finalize()
  ELSE
     CALL COCG_finalize()
  END IF
  !
  ! Output to a file
  !
  IF(ndim == nl) THEN
     CALL output_result_debug()
  ELSE
     CALL output_result()
  END IF
  !
  DEALLOCATE(v12, v2, r_l, x, z)
  IF(lBiCG) DEALLOCATE(v14, v4)
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Done  #####"
  WRITE(*,*)
  !
END PROGRAM shiftk
