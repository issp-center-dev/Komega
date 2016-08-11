MODULE shiftk_vals
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & ndim,    & ! Size of Hilvert space
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
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & z(:)         ! (nomega): Frequency
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & ham(:,:), &
  & rhs(:), &
  & v12(:), v2(:), & ! (ndim): Working vector
  & r_l(:), & ! (ndim) : Projeccted residual vector 
  & x(:,:) ! (ndim,nomega) : Projected result 
  !
  ! Variables for Restart
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & alpha(:), beta(:) ! (iter_old) 
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & r_l_save(:,:) ! (ndim,iter_old) Projected residual vectors
  !
  CHARACTER(256),SAVE :: &
  & inham, &
  & invec, &
  & calctype
  !
  LOGICAL,SAVE :: &
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
!
!
SUBROUTINE input_hamiltonian()
  !
  USE shiftk_vals, ONLY : ndim, ham, inham
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ndim2, nelem, ielem, ii, jj
  REAL(8) :: ham_r, ham_i
  CHARACTER(100) :: ctmp
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Hamiltonian  ##########"
  WRITE(*,*)
  !
  OPEN(fi, file = TRIM(inham))
  READ(fi, *) ctmp
  READ(fi,*) ndim, ndim2, nelem
  WRITE(*,*) "          Dim. of Hamiltonian : ", ndim, ndim2
  WRITE(*,*) "  Num. of Non-Zero Components : ", nelem
  !
  IF(ndim2 /= ndim) THEN
     WRITE(*,*) "ERROR ! Hamiltonian is not square."
     STOP
  END IF
  !
  ALLOCATE(ham(ndim, ndim))
  !
  ham(1:ndim, 1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
  DO ielem = 1, nelem
     READ(fi,*) ii, jj, ham_r, ham_i
     ham(ii,jj) = CMPLX(ham_r, ham_i, KIND(0d0))
  END DO
  !
  CLOSE(fi)
  !
END SUBROUTINE input_hamiltonian
!
!
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
  USE shiftk_vals, ONLY : iter_old, alpha, beta, z_seed, ndim, r_l_save
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, iter, idim
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
  ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(ndim, iter_old))
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
     DO idim = 1, ndim
        READ(fi,*) r_l_save_r, r_l_save_i
        r_l_save(idim,iter) = CMPLX(r_l_save_r, r_l_save_i, KIND(0d0))
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
  USE shiftk_vals, ONLY : ndim, v2, v12
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
  CLOSE(fi)
  !
END SUBROUTINE input_restart_vector
!
! Input restart variables from file
!
SUBROUTINE output_restart_parameter()
  !
  USE shiftk_vals, ONLY : iter_old, alpha, beta, z_seed, r_l_save, ndim
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, iter, idim
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
     DO idim = 1, ndim
        WRITE(fo,'(2e25.16)') r_l_save(idim,iter)
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
  USE shiftk_vals, ONLY : ndim, v2, v12
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
  CLOSE(fo)
  !
END SUBROUTINE output_restart_vector
!
! Check Result
!
SUBROUTINE output_result()
  !
  USE mathlib, ONLY : dgemv
  USE shiftk_vals, ONLY : v2, ndim, x, rhs, z, nomega, ham
  !
  IMPLICIT NONE
  !
  INTEGER :: iz
  CHARACTER(100) :: cndim, form
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
     v2(1:ndim) = z(iz) * x(1:ndim,iz) - rhs(1:ndim)
     CALL zgemv("N", ndim, ndim, CMPLX(-1d0, 0d0, KIND(0d0)), Ham, ndim, x(1:ndim,iz), 1, CMPLX(1d0, 0d0, KIND(0d0)), v2, 1)
     !
     !write(*,form) v2(1:ndim)
     write(*,'(a,i5,a,2e13.5,a,e13.5,a,2e13.5)') "DEBUG (", iz, "), omega = ", z(iz), &
     &            ", Res. = ", dble(dot_product(v2,v2)), ", <rhs|x> = ", dot_product(rhs,x(1:ndim,iz))
     !
  END DO
  !
END SUBROUTINE output_result
!
END MODULE shiftk_routines
!
!
!
PROGRAM shiftk
  !
  USE shifted_cocg, ONLY : COCG_init, COCG_restart, COCG_update, &
  &                        COCG_getcoef, COCG_getvec, COCG_finalize
  USE shiftk_routines, ONLY : input_filename, input_hamiltonian, input_rhs_vector, &
  &                           input_parameter, input_restart_parameter, input_restart_vector, &
  &                           output_result, output_restart_parameter, output_restart_vector
  USE shiftk_vals, ONLY : alpha, beta, calctype, ndim, nomega, maxloops, iter_old, ham, &
  &                       outrestart, rhs, v12, v2, r_l, r_l_save, threshold, x, z, z_seed
  USE mathlib, ONLY : zgemv
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
  CALL input_parameter()
  !
  ALLOCATE(v12(ndim), v2(ndim), r_l(ndim), x(ndim,nomega))
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
        CALL COCG_restart(ndim, ndim, nomega, x, z, maxloops, threshold, status, &
        &                 iter_old, v2, v12, alpha, beta, z_seed, r_l_save)
     ELSE
        CALL COCG_restart(ndim, ndim, nomega, x, z, 0,        threshold, status, &
        &                 iter_old, v2, v12, alpha, beta, z_seed, r_l_save)
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
        CALL COCG_init(ndim, ndim, nomega, x, z, maxloops, threshold, status)
     ELSE
        CALL COCG_init(ndim, ndim, nomega, x, z, 0,        threshold, status)
     END IF
     !
  ELSE
     !
     WRITE(*,*) "ERROR ! calctype = ", TRIM(calctype)
     STOP
     !
  END IF
  !
  ! COCG Loop
  !
  WRITE(*,*)
  WRITE(*,*) "#####  CG Iteration  #####"
  WRITE(*,*)
  !
  DO iter = 1, maxloops
     !
     ! Projection of Residual vector into the space
     ! spaned by left vectors
     !
     r_l(1:ndim) = v2(1:ndim)
     !
     ! Matrix-vector product
     !
     CALL zgemv("N", ndim, ndim, CMPLX(1d0, 0d0, KIND(0d0)), Ham, ndim, v2, 1, CMPLX(0d0, 0d0, KIND(0d0)), v12, 1)
     !
     ! Update result x with COCG
     !
     CALL COCG_update(v12, v2, x, r_l, status)
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
     ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(ndim, iter_old))
     !
     CALL COCG_getcoef(alpha, beta, z_seed, r_l_save)
     CALL COCG_getvec(v12)
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
  CALL COCG_finalize()
  !
  ! Output to a file
  !
  CALL output_result()
  !
  DEALLOCATE(v12, v2, r_l, x, z)
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Done  #####"
  WRITE(*,*)
  !
END PROGRAM shiftk
