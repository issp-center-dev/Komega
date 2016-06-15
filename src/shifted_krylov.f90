MODULE shifted_krylov_vals
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
  & pi_save(:), &
  & p(:,:), &
  & v3(:), &
  & alpha_save(:), &
  & beta_save(:), &
  & r_l_save(:,:)
  !
END MODULE shifted_krylov_vals
!
! Routines for real-valiable CG
!
MODULE shifted_krylov
  !
CONTAINS
!
! Shifted Part
!
SUBROUTINE CG_R_shiftedeqn(r_l, x)
  !
  USE shifted_krylov_vals, ONLY : alpha, beta, itermax, nl, nz, &
  &                               p, pi, pi_old, pi_save, z, zseed
  USE mathlib, ONLY : daxpy
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: r_l(nl)
  REAL(8),INTENT(INOUT) :: x(nl,nz)
  !
  INTEGER :: iz
  REAL(8) :: pi_new
  !
  DO iz = 1, nz
     !
     pi_new = (1d0 + alpha * (z(iz) - zseed)) * pi(iz) &
     &      - alpha * beta / alpha_old * (pi_old(iz) - pi(iz))
     p(1:nl,iz) = r_l(1:nl) / pi(iz) &
     &          + (pi_old(iz) / pi(iz))**2 * beta * p(1:nl,iz)
     CALL daxpy(ndim,pi(iz)/ pi_old(iz) * alpha,p(1:nl,iz),1,x(1:nl,iz),1)
     pi_old(iz) = pi(iz)
     pi(iz) = pi_new
     !
     IF(itermax > 0) pi_save(iz,iter) = pi_new
     !
  END DO
  !
END SUBROUTINE CG_R_shiftedeqn
!
! Seed Switching
!
SUBROUTINE CG_R_seed_switch(v2)
  !
  USE mathlib, ONLY : dscal
  USE shifted_krylov_vals, ONLY : alpha, alpha_save, beta_save, &
  &                               iter, itermax, ndim, nz, pi, pi_old, rho, &
  &                               v3, zseed
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(INOUT) :: v3(ndim)
  !
  INTEGER :: izseed
  REAL(8) :: scale
  !
  izseed = MINLOC(pi(1:nz), 1)
  !
  IF(ABS(zseed - z(izseed)) > 1d-12) THEN
     !
     zseed = z(izseed)
     !
     alpha = alpha * pi_old(izseed) / pi(izseed)
     rho = rho / pi_old(izseed)
     !
     scale = 1d0 / pi(izseed)
     CALL dscal(ndim,scale,v2,1)
     CALL dscal(nz,scale,pi,1)
     scale = 1d0 / pi_old(izseed)
     CALL dscal(ndim,scale,v3,1)
     CALL dscal(nz,scale,pi_old,1)
     !
     ! For restarting
     !
     IF(itermax > 0) THEN
        !
        do jter = 1, iter
           alpha_save(jter) = alpha_save(jter) &
           &                * pi_save(izseed, jter - 1) / pi_save(izseed,jter) 
           beta_save(jter) = beta_save(jter) &
           &               * (pi_save(izseed, jter - 2) / pi_save(izseed,jter - 1))**2 
           scale = 1d0 / pi_save(izseed, jter)
           CALL dscal(nz,scale,pi_save(1:nz,jter),1)
        end do
        !
     END IF
     !
  END IF
  !
END SUBROUTINE CG_R_seed_switch
!
! Allocate & initialize variables
!
SUBROUTINE CG_R_init(ndim0, nl0, nz0, v2, x, z0, itermax0, threshold0, lconverged)
  !
  USE shifted_krylov, ONLY : alpha, alpha_save, beta, beta_save, iter, itermax, &
  &                          ndim, nl, nz, p, pi, pi_old, pi_save, rho, r_l_save, v3, z, zseed 
  USE mathlib, ONLY : dcopy
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ndim0, nl0, nz0, itermax0
  REAL(8),INTENT(IN) :: z0(nz0), threshold0
  REAL(8),INTENT(OUT) :: x(nl0,nz0)
  !
  ndim = ndim0
  nl = nl0
  nz = nz0
  itermax = itermax0
  threshold = threshold0
  !
  ALLOCATE(z(nz), v3(ndim), pi(nz), pi_old(nz), p(nl,nz))
  CALL dcopy(nz,z0,1,z,1)
  v3(1:ndim) = 0d0
  p(1:nl,1:nz) = 0d0
  x(1:nl,1:nz) = 0d0
  pi(1:nz) = 1d0
  pi_old(1:nz) = 1d0
  rho = 1d0
  alpha = 1d0
  beta = 0d0
  zseed = 0d0
  iter = 0
  !
  IF(itermax > 0) THEN
     ALLOCATE(alpha_save(itermax), beta_save(itermax), &
     &        r_l_save(nl,itermax), pi_save(nz,-1:itermax))
     pi_save(1:nz,-1:0) = 1d0
  END IF
  !
  ! Convergence check
  !
  res = ddot(ndim,v2,1,v2,1)
  !
  IF(res < threshold) THEN
     lconverged = 1
  ELSE
     lconverged = 0
  END IF
  !
END SUBROUTINE CG_R_init
!
! Restart by input
!
SUBROUTINE CG_R_restart(ndim0, nl0, nz0, v2, x, z0, itermax0, threshold0, lconverged, &
&                       iter_old, v12, alpha_save0, beta_save0, zseed0, r_l_save0)
  !
  USE shifted_krylov_vals, ONLY : alpha, alpha_old, alpha_save, beta, beta_save, iter, iter_max, &
  &                               ndim, nl, r_l_save, threshold, v3, zseed
  USE mathlib, ONLY : dcopy, daxpy
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ndim0, nl0, nz0, itermax0
  REAL(8),INTENT(IN) :: z0(nz0), threshold0
  REAL(8),INTENT(OUT) :: x(nl0,nz0)
  INTEGER,INTENT(OUT) :: lconverged
  !
  ! For Restarting
  !
  INTEGER(8),INTENT(IN) :: iter_old
  REAL(8),INTENT(IN) :: &
  & r(ndim), r_old(ndim), alpha_save0(iter_old), beta_save0(iter_old), &
  & zseed0, r_l_save0(nl,iter_old)
  !
  CALL CG_R_init(ndim0, nl0, nz0, v2, x, z0, itermax0, threshold0, lconverged)
  zseed = zseed0
  !
  DO iter = 1, iter_old
     !
     beta = beta_save0(iter)
     alpha_old = alpha
     alpha = alpha_save0(iter)
     !
     ! For restarting
     !
     IF(itermax > 0) THEN
        alpha_save(iter) = alpha
        beta_save(iter) = beta
        CALL dcopy(nl,r_l_save0(1:nl,iter),1,r_l_save(1:nl,iter),1)
     END IF
     !
     ! Shifted equation
     !
     CALL CG_R_shiftedeqn(r_l_save0(1:nl,iter), x)
     !
  END DO
  !
  ! Rewind
  !
  iter = iter_old 
  !
  CALL dcopy(ndim,v12,1,v3,1)
  !
  ! Seed Switching
  !
  CALL CG_R_seed_switch(v2)
  !
  ! Convergence check
  !
  res = ddot(ndim,v2,1,v2,1)
  !
  IF(res < threshold) THEN
     lconverged = iter
  ELSE IF(iter == itermax)
     lconverged = - iter
  ELSE
     lconverged = 0
  END IF
  !
END SUBROUTINE CG_R_restart
!
! Update x, p, r
!
SUBROUTINE CG_R_update(v12, v2, x, r_l, lconverged)
  !
  USE mathlib, ONLY : daxpy, ddot, dcopy, dscale
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(INOUT) :: v12(ndim), v2(ndim), x(nl,nz)
  REAL(8),INTENT(IN) :: r_l(nl)
  INTEGER,INTENT(OUT) :: lconverged
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
     CALL dcopy(nl,r_l,1,r_l_save(1:nl,iter),1)
  END IF
  !
  ! Shifted equation
  !
  CALL CG_R_shiftedeqn(r_l, x)
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
  CALL CG_R_seed_switch(v2)
  !
  ! Convergence check
  !
  res = ddot(ndim,v2,1,v12,1)
  IF(res < threshold) THEN
     lconverged = iter
  ELSE IF(iter == itermax)
     lconverged = -iter
  ELSE
     lconverged = 0
  END IF
  !
END SUBROUTINE CG_R_update
!
! Return saved alpha, beta, r_l
!
SUBROUTINE CG_R_getcoef(alpha_save0, beta_save0, zseed0, r_l_save0)
  !
  USE shifted_krylov, ONLY : alpha_save, beta_save, r_l_save, iter, nl, zseed
  USE mathlib, ONLY : dcopy
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(OUT) :: alpha0(iter), beta(iter), zseed0, &
  &                      r_l_save0(nl,iter)
  !
  zseed0 = zseed
  CALL dcopy(iter,alpha_save,1,alpha_save0,1)
  CALL dcopy(iter,beta_save,1,beta_save0,1)
  CALL dcopy(nl*iter,r_l_save,1,r_l_save0,1)
  !
END SUBROUTINE CG_R_getcoef
!
! Return r_old
!
SUBROUTINE CG_R_getvec(r_old)
  !
  USE shifted_krylof, ONLY : ndim, v3
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
! Deallocate private arrays
!
SUBROUTINE CG_R_finalize()
  !
  USE shifted_krylov, ONLY : alpha_save, beta_save, r_l_save, pi_save
  !
  IMPLICIT NONE
  !
  DEALLOCATE(z, v3, pi, pi_old, p)
  !
  IF(itermax > 0) THEN
     DEALLOCATE(alpha_save, beta_save, r_l_save, pi_save)
  END IF
  !
END SUBROUTINE CG_R_finalize
!
END MODULE syfted_krylov
