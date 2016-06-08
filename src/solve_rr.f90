PROGRAM solve_rr
  !
  USE shifted_krilov, ONLY : CG_R_init, CG_R_restart, CG_R_update, &
  &                          CG_R_getcoef, CG_R_getvec, CG_R_finalize
  USE my_module, ONLY : input_size, rhs_vector, input_restart, &
  &                     projection, matrix_product, &
  &                     output_restart, output_result
  !
  IMPLICIT NONE
  !
  INTEGER :: &
  & ndim,    & ! Size of Hilvert space
  & nz,      & ! Number of frequencies
  & nl,      & ! Number of Left vector
  & itermax, & ! Max. number of iteraction
  & itermin, & ! First iteration in this run
  & iter,    & ! Counter for Iteration
  & niter_old  ! Number of iteraction of previous run
  !
  REAL(8) :: &
  & threshold, & ! Convergence Threshold
  !
  REAL(8) :: &
  & v12(:), v2(:), & ! (ndim): Working vector
  & rsmall(:), & ! (nl) : Projeccted residual vector 
  & x(:,:),    & ! (nl,nz) : Projected result 
  & z(:)         ! (nz): Frequency
  !
  ! Variables for Restart
  !
  REAL(8) :: &
  & r(:), r_old(:), & ! (ndim) Residual vector in Current & Previous step
  & zseed, & ! Seed frequency
  & alpha(:), beta(:), & ! (niter_old) 
  & rsmall_save(:,:) ! (nl,niter_old) Projected residual vectors
  !
  LOGICAL :: lconverged ! .TRUE. if converged
  !
  ! Input Size of vectors
  !
  CALL input_size(ndim,nz,nl,itermax)
  !
  ALLOCATE(v12(ndim), v2(ndim), rsmall(nl), x(nl,nz), z(nz))
  !
  ! Generate Right Hand Side Vector
  !
  CALL RANDOM_NUMBER(v2(1:ndim))
  !
  ! Check: Whether the restart file is exist.
  !
  CALL input_restart(niter_old, r, r_old, alpha, beta, zseed, rsmall_save)
  !
  IF(niter_old > 0) THEN
    !
    ! When restarting, counter
    !
    itermin = niter_old + 1
    CALL CG_R_restart(ndim, nl, nz, v2, x, z, itermax, &
    &                 threshold, r, r_old, &
    &                 alpha, beta, zseed, rsmall_save)
    !
    ! These vectors were saved in CG_R routine
    !
    DEALLOCATE(r, r_old, alpha, beta, zseed, rsmall_save)
    !
  ELSE
    itermin = 1
    CALL CG_R_init(ndim, nl, nz, v2, x, z, itermax, threshold)
  END IF
  !
  ! CG_R Loop
  !
  DO iter = itermin, itermax
     !
     ! Projection of Residual vector into the space
     ! spaned by left vectors
     !
     rsmall(1:nl) = v2(1:ndim)
     !
     ! Matrix-vector product
     !
     CALL dgemv("N", ndim, ndim, 1d0, Ham, ndim, v2, 1, 0d0, v12, 1)
     !
     ! Update result x with CG_R
     !
     CALL CG_R_update(v12, v2, x, rsmall, lconverged)
     !
     IF(lconverged == .true.) EXIT
     !
  END DO
  !
  ! Get these vectors for restart in the Next run
  !
  ALLOCATE(r(ndim), r_old(ndim), alpha(iter), beta(iter), rsmall_save(nl,iter))
  !
  CALL CG_R_getcoef(alpha, beta, zseed, rsmall_save)
  CALL CG_R_getvec(r, r_old)
  !
  ! Deallocate all intrinsic vectors
  !
  CALL CG_R_finalize()
  !
  ! Output to a file
  !
  CALL output_restart(niter_old, r, r_old, alpha, beta, zseed, rsmall_save)
  !
  CALL output_result(z,x)
  !
  DEALLOCATE(r, r_old, alpha, beta, zseed, rsmall_save)
  DEALLOCATE(v12, v2, rsmall, x, z)
  !
END PROGRAM solve_rr
