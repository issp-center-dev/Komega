!
! Routines for real-valiable CG
!
MODULE shifted_bicg
  !
CONTAINS
!
! Shifted Part
!
SUBROUTINE BiCG_shiftedeqn(r_l, x)
  !
  USE shifted_krylov_parameter, ONLY : iter, itermax, nl, nz
  USE shifted_krylov_vals_c, ONLY : alpha, alpha_old, beta, pi, pi_old, pi_save, z, z_seed
  USE shifted_krylov_vecs_c, ONLY : p
  USE shifted_krylov_math, ONLY : zaxpy
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(IN) :: r_l(nl)
  COMPLEX(8),INTENT(INOUT) :: x(nl,nz)
  !
  INTEGER :: iz
  COMPLEX(8) :: pi_new
  !
  DO iz = 1, nz
     !
     pi_new = (1d0 + alpha * (z(iz) - z_seed)) * pi(iz) &
     &      - alpha * beta / alpha_old * (pi_old(iz) - pi(iz))
     p(1:nl,iz) = r_l(1:nl) / pi(iz) &
     &          + (pi_old(iz) / pi(iz))**2 * beta * p(1:nl,iz)
     CALL zaxpy(nl, pi(iz)/ pi_new * alpha, p(1:nl,iz), 1, x(1:nl,iz), 1)
     pi_old(iz) = pi(iz)
     pi(iz) = pi_new
     !
     IF(itermax > 0) pi_save(iz,iter) = pi_new
     !
  END DO
  !
END SUBROUTINE BiCG_shiftedeqn
!
! Seed Switching
!
SUBROUTINE BiCG_seed_switch(v2, v4, status)
  !
  USE shifted_krylov_parameter, ONLY : iter, itermax, ndim, nz, nl, threshold
  USE shifted_krylov_vals_c, ONLY : alpha, alpha_save, beta_save, pi, pi_old, &
  &                               pi_save, rho, z, z_seed
  USE shifted_krylov_vecs_c, ONLY : v3, v5, r_l_save
  USE shifted_krylov_math, ONLY : dscal, zscal
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(INOUT) :: v2(ndim), v4(ndim)
  INTEGER,INTENT(OUT) :: status
  !
  INTEGER :: iz_seed, jter
  COMPLEX(8) :: scale
  !
  iz_seed = MINLOC(ABS(pi(1:nz)), 1)
  status = iz_seed
  !
  IF(ABS(z_seed - z(iz_seed)) > 1d-12) THEN
     !
     z_seed = z(iz_seed)
     !
     IF(ABS(pi(iz_seed)) < threshold) THEN
        status = - iz_seed
     ELSE
        status = iz_seed
     END IF
     !
     alpha = alpha * pi_old(iz_seed) / pi(iz_seed)
     rho = rho / pi_old(iz_seed)
     !
     scale = 1d0 / pi(iz_seed)
     CALL zscal(ndim, scale, v2, 1)
     scale = 1d0 / CONJG(pi(iz_seed))
     CALL zscal(ndim, scale, v4, 1)
     !
     scale = 1d0 / pi(iz_seed)
     CALL zscal(nz,scale,pi,1)
     !
     scale = 1d0 / pi_old(iz_seed)
     CALL zscal(ndim, scale, v3, 1)
     scale = 1d0 / CONJG(pi_old(iz_seed))
     CALL zscal(ndim, scale, v5, 1)
     !
     scale = 1d0 / pi_old(iz_seed)
     CALL zscal(nz,scale,pi_old,1)
     !
     ! For restarting
     !
     IF(itermax > 0) THEN
        !
        DO jter = 1, iter
           !
           alpha_save(jter) = alpha_save(jter) &
           &                * pi_save(iz_seed, jter - 1) / pi_save(iz_seed,jter) 
           beta_save(jter) = beta_save(jter) &
           &               * (pi_save(iz_seed, jter - 2) / pi_save(iz_seed,jter - 1))**2 
           !
           scale = 1d0 / pi_save(iz_seed, jter - 1)
           CALL zscal(nl, scale, r_l_save(1:nl,jter), 1)
           !
        END DO
        !
        DO jter = 1, iter
           scale = 1d0 / pi_save(iz_seed, jter)
           CALL zscal(nz,scale,pi_save(1:nz,jter),1)
        END DO
        !
     END IF
     !
  END IF
  !
END SUBROUTINE BiCG_seed_switch
!
! Allocate & initialize variables
!
SUBROUTINE BiCG_init(ndim0, nl0, nz0, x, z0, itermax0, threshold0, status)
  !
  USE shifted_krylov_parameter, ONLY : iter, itermax, ndim, nl, nz, threshold
  USE shifted_krylov_vals_c, ONLY : alpha, alpha_save, beta, beta_save, pi, &
  &                               pi_old, pi_save, rho, z, z_seed 
  USE shifted_krylov_vecs_c, ONLY : p, r_l_save, v3, v5
  USE shifted_krylov_math, ONLY : zcopy
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ndim0, nl0, nz0, itermax0
  REAL(8),INTENT(IN) :: threshold0
  COMPLEX(8),INTENT(IN) :: z0(nz0)
  COMPLEX(8),INTENT(OUT) :: x(nl0,nz0)
  INTEGER,INTENT(OUT) :: status(3)
  !
  status(1:3) = 0
  !
  ndim = ndim0
  nl = nl0
  nz = nz0
  itermax = itermax0
  threshold = threshold0
  !
  ALLOCATE(z(nz), v3(ndim), v5(ndim), pi(nz), pi_old(nz), p(nl,nz))
  CALL zcopy(nz,z0,1,z,1)
  v3(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
  v5(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
  p(1:nl,1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  x(1:nl,1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  pi(1:nz) = CMPLX(1d0, 0d0, KIND(0d0))
  pi_old(1:nz) = CMPLX(1d0, 0d0, KIND(0d0))
  rho = CMPLX(1d0, 0d0, KIND(0d0))
  alpha = CMPLX(1d0, 0d0, KIND(0d0))
  beta = CMPLX(0d0, 0d0, KIND(0d0))
  z_seed = z(1)
  iter = 0
  !
  IF(itermax > 0) THEN
     ALLOCATE(alpha_save(itermax), beta_save(itermax), &
     &        r_l_save(nl,itermax), pi_save(nz,-1:itermax))
     pi_save(1:nz,-1:0) = CMPLX(1d0, 0d0, KIND(0d0))
  END IF
  !
END SUBROUTINE BiCG_init
!
! Restart by input
!
SUBROUTINE BiCG_restart(ndim0, nl0, nz0, x, z0, itermax0, threshold0, status, &
&                       iter_old, v2, v12, v4, v14, alpha_save0, beta_save0, z_seed0, r_l_save0)
  !
  USE shifted_krylov_parameter, ONLY : iter, itermax, ndim, nl, threshold
  USE shifted_krylov_vals_c, ONLY : alpha, alpha_old, alpha_save, beta, beta_save, z_seed
  USE shifted_krylov_vecs_c, ONLY : r_l_save, v3, v5
  USE shifted_krylov_math, ONLY : zcopy, zdotc
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ndim0, nl0, nz0, itermax0
  REAL(8),INTENT(IN) :: threshold0
  COMPLEX(8),INTENT(IN) :: z0(nz0)
  COMPLEX(8),INTENT(OUT) :: x(nl0,nz0)
  INTEGER,INTENT(OUT) :: status(3)
  !
  ! For Restarting
  !
  INTEGER,INTENT(IN) :: iter_old
  COMPLEX(8),INTENT(IN) :: &
  & alpha_save0(iter_old), beta_save0(iter_old), z_seed0
  COMPLEX(8),INTENT(IN) :: r_l_save0(nl0,iter_old)
  COMPLEX(8),INTENT(INOUT) :: v2(ndim), v12(ndim)
  COMPLEX(8),INTENT(INOUT) :: v4(ndim), v14(ndim)
  !
  CALL BiCG_init(ndim0, nl0, nz0, x, z0, itermax0, threshold0, status)
  z_seed = z_seed0
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
        CALL zcopy(nl,r_l_save0(1:nl,iter),1,r_l_save(1:nl,iter),1)
     END IF
     !
     ! Shifted equation
     !
     CALL BiCG_shiftedeqn(r_l_save0(1:nl,iter), x)
     !
  END DO
  !
  ! Rewind
  !
  iter = iter_old 
  !
  CALL zcopy(ndim,v12,1,v3,1)
  CALL zcopy(ndim,v14,1,v5,1)
  !
  ! Seed Switching
  !
  CALL BiCG_seed_switch(v2,v4,status(3))
  !
  ! Convergence check
  !
  v12(1) = CMPLX(MAXVAL(ABS(v2(1:ndim))), 0d0, KIND(0d0))
  !
  IF(DBLE(v12(1)) < threshold) THEN
     status(1) = iter
  ELSE IF(iter == itermax) THEN
     status(1) = - iter
  ELSE
     status(1) = 0
  END IF
  !
END SUBROUTINE BiCG_restart
!
! Update x, p, r
!
SUBROUTINE BiCG_update(v12, v2, v14, v4, x, r_l, status)
  !
  USE shifted_krylov_parameter, ONLY : iter, itermax, ndim, nl, nz, threshold
  USE shifted_krylov_vals_c, ONLY : alpha, alpha_old, alpha_save, &
  &                               beta, beta_save, rho, z_seed
  USE shifted_krylov_vecs_c, ONLY : r_l_save, v3, v5
  USE shifted_krylov_math, ONLY : zdotc, zcopy
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(INOUT) :: v12(ndim), v2(ndim), v14(ndim), v4(ndim), x(nl,nz)
  COMPLEX(8),INTENT(IN) :: r_l(nl)
  INTEGER,INTENT(INOUT) :: status(3)
  !
  COMPLEX(8) :: rho_old, alpha_denom
  !
  iter = iter + 1
  !status(1:3) = 0
  !
  rho_old = rho
  rho = zdotc(ndim,v4,1,v2,1)
  IF(iter == 1) THEN
     beta = CMPLX(0d0, 0d0, KIND(0d0))
  ELSE
     beta = rho / rho_old
  END IF
  v12(1:ndim) = z_seed * v2(1:ndim) - v12(1:ndim)
  v14(1:ndim) = CONJG(z_seed) * v4(1:ndim) - v14(1:ndim)
  alpha_old = alpha
  alpha_denom = zdotc(ndim,v4,1,v12,1) - beta * rho / alpha
  !
  IF(ABS(alpha_denom) < threshold) THEN
     status(2) = 1
     alpha = 1d0
  ELSE IF(ABS(rho) < threshold) THEN
     status(2) = 2
     alpha = 1d0
     beta = 0d0
  ELSE
     status(2) = 0
     alpha = rho / alpha_denom
  END IF
  !
  ! For restarting
  !
  IF(itermax > 0) THEN
     alpha_save(iter) = alpha
     beta_save(iter) = beta
     CALL zcopy(nl,r_l,1,r_l_save(1:nl,iter),1)
  END IF
  !
  ! Shifted equation
  !
  CALL BiCG_shiftedeqn(r_l, x)
  !
  ! Update residual
  !
  v12(1:ndim) = (1d0 + alpha * beta / alpha_old) * v2(1:ndim) &
  &           - alpha * v12(1:ndim) &
  &           - alpha * beta / alpha_old * v3(1:ndim)
  CALL zcopy(ndim,v2,1,v3,1)
  CALL zcopy(ndim,v12,1,v2,1)
  v14(1:ndim) = (1d0 + CONJG(alpha * beta / alpha_old)) * v4(1:ndim) &
  &           - CONJG(alpha) * v14(1:ndim) &
  &           - CONJG(alpha * beta / alpha_old) * v5(1:ndim)
  CALL zcopy(ndim,v4,1,v5,1)
  CALL zcopy(ndim,v14,1,v4,1)
  !
  ! Seed Switching
  !
  !CALL BiCG_seed_switch(v2,v4,status(3))
  !
  ! Convergence check
  !
  rho_old = CMPLX(MAXVAL(ABS(v2(1:ndim))), 0d0, KIND(0d0))
  !
  IF(DBLE(rho_old) < threshold) THEN
     status(1) = iter
  ELSE IF(iter == itermax) THEN
     status(1) = -iter
  ELSE
     v12(1) = rho_old
     status(1) = 0
  END IF
  !
END SUBROUTINE BiCG_update
!
! Return saved alpha, beta, r_l
!
SUBROUTINE BiCG_getcoef(alpha_save0, beta_save0, z_seed0, r_l_save0)
  !
  USE shifted_krylov_parameter, ONLY : iter, nl
  USE shifted_krylov_vals_c, ONLY : alpha_save, beta_save, z_seed
  USE shifted_krylov_vecs_c, ONLY : r_l_save
  USE shifted_krylov_math, ONLY : dcopy, zcopy
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(OUT) :: alpha_save0(iter), beta_save0(iter), z_seed0
  COMPLEX(8),INTENT(OUT) :: r_l_save0(nl,iter)
  !
  z_seed0 = z_seed
  CALL zcopy(iter,alpha_save,1,alpha_save0,1)
  CALL zcopy(iter,beta_save,1,beta_save0,1)
  CALL zcopy(nl*iter,r_l_save,1,r_l_save0,1)
  !
END SUBROUTINE BiCG_getcoef
!
! Return r_old
!
SUBROUTINE BiCG_getvec(r_old, r_tilde_old)
  !
  USE shifted_krylov_parameter, ONLY : ndim
  USE shifted_krylov_vecs_c, ONLY : v3, v5
  USE shifted_krylov_math, ONLY : zcopy
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(OUT) :: r_old(ndim), r_tilde_old(ndim)
  !
  CALL zcopy(ndim,v3,1,r_old,1)
  CALL zcopy(ndim,v5,1,r_tilde_old,1)
  !
END SUBROUTINE BiCG_getvec
!
! Deallocate private arrays
!
SUBROUTINE BiCG_finalize()
  !
  USE shifted_krylov_parameter, ONLY : itermax
  USE shifted_krylov_vals_c, ONLY : alpha_save, beta_save, &
  &                                 pi, pi_old, pi_save, z
  USE shifted_krylov_vecs_c, ONLY : p, r_l_save, v3, v5
  !
  IMPLICIT NONE
  !
  DEALLOCATE(z, v3, v5, pi, pi_old, p)
  !
  IF(itermax > 0) THEN
     DEALLOCATE(alpha_save, beta_save, r_l_save, pi_save)
  END IF
  !
END SUBROUTINE BiCG_finalize
!
END MODULE shifted_bicg
