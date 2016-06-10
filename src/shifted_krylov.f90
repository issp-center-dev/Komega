MODULE syfted_krylov
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & ndim, &
  & nl,   &
  & nz,   &
  & itermax, &
  & iter
  !
  REAL(8),SAVE :: &
  & threshold, &
  & zseed, &
  & rho, &
  & alpha, &
  & beta
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & z(:), &
  & v3(:), &
  & pi(:), &
  & pi_old(:), &
  & p(:,:), &
  & v3(:), &
  & alpha_save(:), &
  & beta_save(:), &
  & r_small_save(:,:)
  !
CONTAINS
!
SUBROUTINE CG_R_restart(ndim0, nl0, nz0, v2, x, z0, itermax0, &
&                       threshold0, r, r_old, &
&                       alpha_save0, beta_save0, zseed0, r_small_save0)
  !
  USE mathlib, ONLY : dcopy
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ndim0, nl0, nz0
  REAL(8),INTENT(IN) :: z0(nz0)
  REAL(8),INTENT(OUT) :: x(nl0,nz0)
  INTEGER,INTENT(IN) :: itermax0, threshold0
  !
  ! For Restarting
  !
  REAL(8),INTENT(IN) :: &
  & r(ndim), r_old(ndim), alpha_save0(iter_old), beta_save0(iter_old), &
  & zseed0, r_small_save0(nl,iter_old)
  !
  CALL CG_R_init(ndim0, nl0, nz0, x, z0, itermax0, threshold0)


END SUBROUTINE CG_R_restart
!
SUBROUTINE CG_R_init(ndim0, nl0, nz0, x, z0, itermax0, threshold0)
  !
  USE mathlib, ONLY : dcopy
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ndim0, nl0, nz0
  REAL(8),INTENT(IN) :: z0(nz0)
  REAL(8),INTENT(OUT) :: x(nl0,nz0)
  INTEGER,INTENT(IN) :: itermax0, threshold0
  !
  ndim = ndim0
  nl = nl0
  nz = nz0
  itermax = itermax0
  threshold = threshold0
  !
  ALLOCATE(z(nz), v3(ndim), pi(nz), piold(nz), p(nl,nz))
  CALL dcopy(nz,z0,1,z,1)
  v3(1:ndim) = 0d0
  p(1:nl,1:nz) = 0d0
  x(1:nl,1:nz) = 0d0
  pi(1:nz) = 1d0
  piold(1:nz) = 1d0
  rho = 1d0
  alpha = 1d0
  beta = 0d0
  zseed = 0d0
  iter = 0
  !
  IF(itermax > 0) THEN
     ALLOCATE(alpha_save(itermax), beta_save(itermax), &
     &        r_small_save(1:nl,itermax))
  END IF
  !
END SUBROUTINE CG_R_init
!
SUBROUTINE CG_R_update(v12, v2, x, r_small, lconverged)
  !
  USE mathlib, ONLY : daxpy, ddot, dcopy, dscale
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(INOUT) :: v12(ndim), v2(ndim), x(nl,nz)
  REAL(8),INTENT(IN) :: r_small(nl)
  LOGICAL,INTENT(OUT) :: lconverged
  !
  INTEGER :: iz, izseed
  REAL(8) :: rho_old, prod, alpha_old, scale, res
  !
  iter = iter + 1
  !
  rho_old = rho
  rho = ddot(ndim,v2,1,v2,1)
  beta = rho / rho_old
  v12(1:ndim) = zseed * v2(1:ndim) - v12(1:ndim)
  alpha_old = alpha
  prod = ddot(ndim,v2,1,v12,1)
  alpha = rho / (prod - beta * rho / alpha)
  !
  ! For restarting
  !
  IF(itermax > 0) THEN
     alpha_save(iter) = alpha
     beta_save(iter) = beta
     CALL dcopy(nl,r_small,1,r_small_save(1:nl,iter),1)
  END IF
  !
  ! Shifted equation
  !
  DO iz = 1, nz
     pi_new = (1d0 + alpha * (z(iz) - zseed)) * pi(iz) &
     &      - alpha * beta / alpha_old * (pi_old(iz) - pi(iz))
     p(1:nl,iz) = r_small_(1:nl) / pi(iz) &
     &          + (pi_old(iz) / pi(iz))**2 * beta * p(1:nl,iz)
     CALL daxpy(ndim,pi(iz)/ pi_old(iz) * alpha,p(1:nl,iz),1,x(1:nl,iz),1)
  END DO
  !
  ! Update residual
  !
  v12(1:ndim) = (1d0 + alpha * beta / alpha_old) * v2(1:ndim) &
  &           - alpha * v12(1:ndim)
  &           - alpha * beta / alpha_old * v3(1:ndim)
  CALL dcopy(ndim,v2,1,v3,1)
  CALL dcopy(ndim,v12,1,v2,1)
  !
  ! Seed Switching
  !
  izseed = MINLOC(pi(1:nz), 1)
  zseed = z(izseed)
  scale = 1d0 / pi(izseed)
  CALL dscal(nz,scale,pi,1)
  CALL dscal(ndim,scale,v2,1)
  scale = 1d0 / pi_old(izseed)
  CALL dscal(nz,scale,pi_old,1)
  CALL dscal(ndim,scale,v3,1)
  !
  ! Convergence check
  !
  res = ddot(ndim,v2,1,v12,1)
  IF(res < threshold) THEN
     lconverged = .TRUE.
  ELSE
     lconverged = .FALSE.
  END IF
  !
END SUBROUTINE CG_R_update
!
SUBROUTINE CG_R_getcoef(alpha_save0, beta_save0, zseed0, r_small_save0)
  !
  USE mathlib, ONLY : dcopy
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(OUT) :: alpha0(iter), beta(iter), zseed0, &
  &                      r_small_save0(nl,iter)
  !
  zseed0 = zseed
  CALL dcopy(iter,alpha_save,1,alpha_save0,1)
  CALL dcopy(iter,beta_save,1,beta_save0,1)
  CALL dcopy(nl*iter,r_small_save,1,r_small_save0,1)
  !
END SUBROUTINE CG_R_getcoef
!
SUBROUTINE CG_R_getvec(r_old)
  !
  USE mathlib, ONLY : dcopy
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(OUT) :: r_old(ndim)
  !
  CALL dcopy(ndim,v3,1,r_old,1)
  !
END SUBROUTINE CG_R_getvec
!
SUBROUTINE CG_R_finalize()
  DEALLOCATE(z, v3, pi, piold, p)
END SUBROUTINE CG_R_finalize

END MODULE syfted_krylov
