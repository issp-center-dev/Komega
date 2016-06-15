MODULE solver_rr_vals
  !
  IMPLICIT NONE
  !
  INTEGER :: &
  & ndim,    & ! Size of Hilvert space
  & nz,      & ! Number of frequencies
  & nl,      & ! Number of Left vector
  & itermax, & ! Max. number of iteraction
  & iter_old   ! Number of iteraction of previous run
  !
  REAL(8) :: &
  & v12(:), v2(:), & ! (ndim): Working vector
  & r_l(:), & ! (nl) : Projeccted residual vector 
  & x(:,:),    & ! (nl,nz) : Projected result 
  & z(:)         ! (nz): Frequency
  !
  ! Variables for Restart
  !
  REAL(8) :: &
  & zseed, & ! Seed frequency
  & alpha(:), beta(:), & ! (iter_old) 
  & r_l_save(:,:) ! (nl,iter_old) Projected residual vectors
  !
END MODULE solver_rr_vals

MODULE solver_rr_routines
  !
  IMPLICIT NONE
  !
CONTAINS
  !
SUBROUTINE input_size()
  !
  USE solver_rr_vals, ONLY : ndim, nz, itermax

  IMPLICIT NONE
  !
  NAMELIST /input/ ndim, nz, itermax
  !
  read(*,input,err=100)
  !
  WRITE(*,*) "   Dimension of vvector : ", ndim
  WRITE(*,*) "  Number of frequencies : ", nz
  WRITE(*,*) "  Number of left vector : ", nl
  WRITE(*,*) "       Max. iteractions : ", itermax
  !
  return
  !
100 write(*,*) "Stop in stdin. reading namelist file"
  !
  stop
  !
END SUBROUTINE input_size
!
! Input restart variables from file
!
SUBROUTINE input_restart()
  !
  USE solver_rr_vals, ONLY : iter_old, v2, v12, alpha, beta, zseed, r_l_save, nl, ndim
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ierr
  !
  open(fi, file = 'restart.dat',status="old", action = 'read',iostat = ierr)
  !
  IF(ierr == 0) THEN
     !
     read(fi,*) iter_old
     !
     ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl,iter_old))
     read(fi,*) zseed, v2(1:ndim), v12(1:ndim), &
     &          alpha(1:iter_old), beta(1:iter_old), r_l_save(1:nl, 1:iter_old)
     !
     close(fi)
     !
  ELSE
     !
     iter_old = 0
     !
  END IF
  !
END SUBROUTINE input_restart
!
! Output variables for restart
!
SUBROUTINE output_restart()
  !
  USE solver_rr_vals, ONLY : iter_old, v2, v12, alpha, beta, zseed, r_l_save, nl, ndim
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20
  !
  open(fo, file = 'restart.dat', action = 'write')
  !
  write(fo,*) iter_old
  write(fo,*) zseed, v2(1:ndim), v12(1:ndim), &
  &          alpha(1:iter_old), beta(1:iter_old), r_l_save(1:nl, 1:iter_old)
  !
END SUBROUTINE input_restart
!
! Check Result
!
SUBROUTINE output_result()
  !
  USE solver_rr_vals, ONLY : iter_old, v2, v12, alpha, beta, zseed, r_l_save, nl, ndim
  !
  IMPLICIT NONE
  !
  DO iz = 1, nz
     !
     v2(1:ndim) = z(iz) * x(1:ndim) - rhs(1:ndim)
     CALL dgemv("N", ndim, ndim, -1d0, Ham, ndim, x(1:ndim,iz), 1, 1d0, v2, 1)
     write(*,*) v2(1:ndim)
     !
  END DO
  !
END SUBROUTINE output_result
!
END MODULE solver_rr_routines
!
!
!
PROGRAM solve_rr
  !
  USE shifted_krilov, ONLY : CG_R_init, CG_R_restart, CG_R_update, &
  &                          CG_R_getcoef, CG_R_getvec, CG_R_finalize
  USE solver_rr_routines, ONLY : input_size, input_restart, &
  &                              output_restart, output_result
  !
  IMPLICIT NONE
  !
  INTEGER :: &
  & itermin, & ! First iteration in this run
  & iter,    & ! Counter for Iteration
  & lconverged
  !
  REAL(8) :: threshold ! Convergence Threshold
  !
  ! Input Size of vectors
  !
  CALL input_size()
  !
  ALLOCATE(v12(ndim), v2(ndim), r_l(nl), x(nl,nz), z(nz), ham(ndim,ndim))
  READ(*,*) z(1:nz)
  CALL RANDOM_NUMBER(ham(1:ndim, 1:ndim))
  CALL RANDOM_NUMBER(rhs(1:ndim))
  !
  ! Check: Whether the restart file is exist.
  !
  CALL input_restart()
  !
  IF(iter_old > 0) THEN
    !
    ! When restarting, counter
    !
    itermin = iter_old + 1
    CALL CG_R_restart(ndim, nl, nz, v2, x, z, itermax, threshold, lconverged, &
    &                 iter_old, v12, alpha_save, beta_save, zseed, r_l_save)
    !
    ! These vectors were saved in CG_R routine
    !
    DEALLOCATE(alpha, beta, r_l_save)
    !
  ELSE
     !
     itermin = 1
     !
     ! Generate Right Hand Side Vector
     !
     v2(1:ndim) = rhs(1:ndim)
     !
     CALL CG_R_init(ndim, nl, nz, v2, x, z, itermax, threshold, lconverged)
     !
  END IF
  !
  ! CG_R Loop
  !
  DO iter = 1, itermax - iter_old
     !
     ! Projection of Residual vector into the space
     ! spaned by left vectors
     !
     r_l(1:nl) = v2(1:ndim)
     !
     ! Matrix-vector product
     !
     CALL dgemv("N", ndim, ndim, 1d0, Ham, ndim, v2, 1, 0d0, v12, 1)
     !
     ! Update result x with CG_R
     !
     CALL CG_R_update(v12, v2, x, r_l, lconverged)
     !
     IF(lconverged /= 0) EXIT
     !
  END DO
  !
  IF(lconverged < 0) THEN
     WRITE(*,*) "  Not Converged in iteration ", -lconverged
  ELSE
     WRITE(*,*) "  Converged in iteration ", lconverged
  END IF
  iter_old = abs(lconverged)
  !
  ! Get these vectors for restart in the Next run
  !
  ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl,iter_old))
  !
  CALL CG_R_getcoef(alpha_save, beta_save, zseed, r_l_save)
  CALL CG_R_getvec(v12)
  !
  ! Deallocate all intrinsic vectors
  !
  CALL CG_R_finalize()
  !
  ! Output to a file
  !
  CALL output_restart()
  !
  CALL output_result()
  !
  DEALLOCATE(alpha, beta, zseed, r_l_save)
  DEALLOCATE(v12, v2, r_l, x, z)
  !
END PROGRAM solve_rr
