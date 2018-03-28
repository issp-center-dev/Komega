!
! ISSP Math Library - A library for solving linear systems in materials science
! Copyright (C) 2016 Mitsuaki Kawamura
! 
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
! 
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
! 
! For more details, See `COPYING.LESSER' in the root directory of this library.
!
!
! Routines for COCG
!
MODULE komega_cocg
  !
  PRIVATE
  !
  PUBLIC komega_COCG_init, komega_COCG_restart, komega_COCG_update, komega_COCG_getcoef, &
  &      komega_COCG_getvec, komega_COCG_finalize, komega_COCG_getresidual
  !
CONTAINS
!
! Shifted Part
!
SUBROUTINE komega_COCG_shiftedeqn(r_l, x)
  !
  USE komega_parameter, ONLY : iter, itermax, nl, nz, lz_conv
  USE komega_vals_c, ONLY : alpha, alpha_old, beta, pi, pi_old, pi_save, z, z_seed
  USE komega_vecs_c, ONLY : p
  USE komega_math, ONLY : zaxpy
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
     IF(lz_conv(iz)) cycle
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
END SUBROUTINE komega_COCG_shiftedeqn
!
! Seed Switching
!
SUBROUTINE komega_COCG_seed_switch(v2,status)
  !
  USE komega_parameter, ONLY : iter, itermax, ndim, nz, nl, iz_seed, almost0, lz_conv
  USE komega_vals_c, ONLY : alpha, alpha_save, beta_save, pi, pi_old, &
  &                               pi_save, rho, z, z_seed
  USE komega_vecs_c, ONLY : v3, r_l_save
  USE komega_math, ONLY : dscal, zscal
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(INOUT) :: v2(ndim)
  INTEGER,INTENT(INOUT) :: status(3)
  !
  INTEGER :: jter
  COMPLEX(8) :: scale
  !
  status(3) = MINLOC(ABS(pi(1:nz)), 1, .NOT. lz_conv(1:nz))
  !
  IF(ABS(pi(status(3))) < almost0) THEN
     status(2) = 3
  END IF
  !
  IF(status(3) /= iz_seed) THEN
     !
     iz_seed = status(3)
     z_seed = z(iz_seed)
     !
     alpha = alpha * pi_old(iz_seed) / pi(iz_seed)
     rho = rho / pi_old(iz_seed)**2
     !
     scale = 1d0 / pi(iz_seed)
     CALL zscal(ndim, scale, v2, 1)
     CALL zscal(nz,scale,pi,1)
     scale = 1d0 / pi_old(iz_seed)
     CALL zscal(ndim, scale, v3, 1)
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
END SUBROUTINE komega_COCG_seed_switch
!
! Allocate & initialize variables
!
SUBROUTINE komega_COCG_init(ndim0, nl0, nz0, x, z0, itermax0, threshold0, comm0) BIND(C)
  !
  USE ISO_C_BINDING
  USE komega_parameter, ONLY : iter, itermax, ndim, nl, nz, &
  &                            threshold, iz_seed, lz_conv, lmpi, comm
  USE komega_vals_c, ONLY : alpha, alpha_save, beta, beta_save, pi, &
  &                               pi_old, pi_save, rho, z, z_seed 
  USE komega_vecs_c, ONLY : p, r_l_save, v3
  USE komega_math, ONLY : dcopy
  !
  IMPLICIT NONE
  !
  INTEGER(C_INT),INTENT(IN) :: ndim0, nl0, nz0, itermax0
  REAL(C_DOUBLE),INTENT(IN) :: threshold0
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(IN) :: z0(nz0)
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(OUT) :: x(nl0,nz0)
  INTEGER(C_INT),INTENT(IN),OPTIONAL :: comm0
  !
  ndim = ndim0
  nl = nl0
  nz = nz0
  itermax = itermax0
  threshold = threshold0
  !
  comm = 0
  IF(PRESENT(comm0)) comm = comm0
  lmpi = .FALSE.
#if defined(__MPI)
  IF(PRESENT(comm0)) lmpi = .TRUE.
#endif
  !
  ALLOCATE(z(nz), v3(ndim), pi(nz), pi_old(nz), p(nl,nz), lz_conv(nz))
  CALL zcopy(nz,z0,1,z,1)
  v3(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
  p(1:nl,1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  x(1:nl,1:nz) = CMPLX(0d0, 0d0, KIND(0d0))
  pi(1:nz) = CMPLX(1d0, 0d0, KIND(0d0))
  pi_old(1:nz) = CMPLX(1d0, 0d0, KIND(0d0))
  rho = CMPLX(1d0, 0d0, KIND(0d0))
  alpha = CMPLX(1d0, 0d0, KIND(0d0))
  beta = CMPLX(0d0, 0d0, KIND(0d0))
  iz_seed = 1
  z_seed = z(iz_seed)
  iter = 0
  lz_conv(1:nz) = .FALSE.
  !
  IF(itermax > 0) THEN
     ALLOCATE(alpha_save(itermax), beta_save(itermax), &
     &        r_l_save(nl,itermax), pi_save(nz,-1:itermax))
     pi_save(1:nz,-1:0) = CMPLX(1d0, 0d0, KIND(0d0))
  END IF
  !
END SUBROUTINE komega_COCG_init
!
! Restart by input
!
SUBROUTINE komega_COCG_restart(ndim0, nl0, nz0, x, z0, itermax0, threshold0, status, &
& iter_old, v2, v12, alpha_save0, beta_save0, z_seed0, r_l_save0, comm0) BIND(C)
  !
  USE ISO_C_BINDING
  USE komega_parameter, ONLY : iter, itermax, ndim, nl, threshold, iz_seed, lz_conv, nz, resnorm
  USE komega_vals_c, ONLY : alpha, alpha_old, alpha_save, beta, beta_save, rho, z_seed, pi
  USE komega_vecs_c, ONLY : r_l_save, v3
  USE komega_math, ONLY : zcopy, zdotuMPI, zdotcMPI
  !
  IMPLICIT NONE
  !
  INTEGER(C_INT),INTENT(IN) :: ndim0, nl0, nz0, itermax0
  REAL(C_DOUBLE),INTENT(IN) :: threshold0
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(IN) :: z0(nz0)
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(OUT) :: x(nl0,nz0)
  INTEGER(C_INT),INTENT(OUT) :: status(3)
  INTEGER(C_INT),INTENT(IN),OPTIONAL :: comm0
  !
  ! For Restarting
  !
  INTEGER(C_INT),INTENT(IN) :: iter_old
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(IN) :: &
  & alpha_save0(iter_old), beta_save0(iter_old), z_seed0
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(IN) :: r_l_save0(nl0,iter_old)
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(INOUT) :: v2(ndim), v12(ndim)
  !
  INTEGER :: iz
  !
  IF(PRESENT(comm0)) THEN
     CALL komega_COCG_init(ndim0, nl0, nz0, x, z0, itermax0, threshold0, comm0)
  ELSE
     CALL komega_COCG_init(ndim0, nl0, nz0, x, z0, itermax0, threshold0)
  END IF
  z_seed = z_seed0
  iz_seed = 0
  !
  status(1:3) = 0
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
     CALL komega_COCG_shiftedeqn(r_l_save0(1:nl,iter), x)
     !
  END DO
  !
  ! Rewind
  !
  iter = iter_old 
  !
  CALL zcopy(ndim,v12,1,v3,1)
  rho = zdotuMPI(ndim,v3,v3)
  !
  ! Seed Switching
  !
  CALL komega_COCG_seed_switch(v2,status)
  !
  ! Convergence check
  !
  v12(1) = CMPLX(SQRT(DBLE(zdotcMPI(ndim,v2,v2))), 0d0, KIND(0d0))
  resnorm = DBLE(v12(1))
  !
  DO iz = 1, nz
     IF(ABS(v12(1)/pi(iz)) < threshold) lz_conv(iz) = .TRUE.
  END DO
  !
  IF(DBLE(v12(1)) < threshold) THEN
     !
     ! Converged
     !
     status(1) = - iter
     status(2) = 0
  ELSE IF(iter == itermax) THEN
     !
     ! NOT Converged in itermax
     !
     status(1) = - iter
     status(2) = 1
  ELSE IF(status(2) == 3) THEN
     !
     ! pi_seed becomes zero
     !
     status(1) = - iter
  ELSE
     !
     ! Continue
     !
     status(1) = iter
     status(2) = 0
  END IF
  !
END SUBROUTINE komega_COCG_restart
!
! Update x, p, r
!
SUBROUTINE komega_COCG_update(v12, v2, x, r_l, status) BIND(C)
  !
  USE ISO_C_BINDING
  USE komega_parameter, ONLY : iter, itermax, ndim, nl, nz, &
  &                            threshold, almost0, lz_conv, resnorm
  USE komega_vals_c, ONLY : alpha, alpha_old, alpha_save, &
  &                         beta, beta_save, rho, z_seed, pi
  USE komega_vecs_c, ONLY : r_l_save, v3
  USE komega_math, ONLY : zdotuMPI, zcopy, zdotcMPI
  !
  IMPLICIT NONE
  !
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(INOUT) :: v12(ndim), v2(ndim), x(nl,nz)
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(IN) :: r_l(nl)
  INTEGER(C_INT),INTENT(INOUT) :: status(3)
  !
  INTEGER :: iz
  COMPLEX(8) :: rho_old, alpha_denom
  !
  iter = iter + 1
  status(1:3) = 0
  !
  rho_old = rho
  rho = zdotuMPI(ndim,v2,v2)
  IF(iter == 1) THEN
     beta = CMPLX(0d0, 0d0, KIND(0d0))
  ELSE
     beta = rho / rho_old
  END IF
  v12(1:ndim) = z_seed * v2(1:ndim) - v12(1:ndim)
  alpha_old = alpha
  alpha_denom = zdotuMPI(ndim,v2,v12) - beta * rho / alpha
  !
  IF(ABS(alpha_denom) < almost0) THEN
     status(2) = 2
  ELSE IF(ABS(rho) < almost0) THEN
     status(2) = 4
  END IF
  alpha = rho / alpha_denom
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
  CALL komega_COCG_shiftedeqn(r_l, x)
  !
  ! Update residual
  !
  v12(1:ndim) = (1d0 + alpha * beta / alpha_old) * v2(1:ndim) &
  &           - alpha * v12(1:ndim) &
  &           - alpha * beta / alpha_old * v3(1:ndim)
  CALL zcopy(ndim,v2,1,v3,1)
  CALL zcopy(ndim,v12,1,v2,1)
  !
  ! Seed Switching
  !
  CALL komega_COCG_seed_switch(v2,status)
  !
  ! Convergence check
  !
  v12(1) = CMPLX(SQRT(DBLE(zdotcMPI(ndim,v2,v2))), 0d0, KIND(0d0))
  resnorm = DBLE(v12(1))
  !
  DO iz = 1, nz
     IF(ABS(v12(1)/pi(iz)) < threshold) lz_conv(iz) = .TRUE.
  END DO
  !
  IF(DBLE(v12(1)) < threshold) THEN
     !
     ! Converged
     !
     status(1) = - iter
     status(2) = 0
  ELSE IF(iter == itermax) THEN
     !
     ! NOT Converged in itermax
     !
     status(1) = - iter
     status(2) = 1
  ELSE IF(status(2) == 2) THEN
     !
     ! alpha becomes infinite
     !
     status(1) = - iter
  ELSE IF(status(2) == 3) THEN
     !
     ! pi_seed becomes zero
     !
     status(1) = - iter
  ELSE IF(status(2) == 4) THEN
     !
     ! rho becomes zero
     !
     status(1) = - iter
  ELSE
     !
     ! Continue
     !
     status(1) = iter
     status(2) = 0
  END IF
  !
END SUBROUTINE komega_COCG_update
!
! Return saved alpha, beta, r_l
!
SUBROUTINE komega_COCG_getcoef(alpha_save0, beta_save0, z_seed0, r_l_save0) BIND(C)
  !
  USE ISO_C_BINDING
  USE komega_parameter, ONLY : iter, nl
  USE komega_vals_c, ONLY : alpha_save, beta_save, z_seed
  USE komega_vecs_c, ONLY : r_l_save
  USE komega_math, ONLY : dcopy, zcopy
  !
  IMPLICIT NONE
  !
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(OUT) :: alpha_save0(iter), beta_save0(iter), z_seed0
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(OUT) :: r_l_save0(nl,iter)
  !
  z_seed0 = z_seed
  CALL zcopy(iter,alpha_save,1,alpha_save0,1)
  CALL zcopy(iter,beta_save,1,beta_save0,1)
  CALL zcopy(nl*iter,r_l_save,1,r_l_save0,1)
  !
END SUBROUTINE komega_COCG_getcoef
!
! Return r_old
!
SUBROUTINE komega_COCG_getvec(r_old) BIND(C)
  !
  USE ISO_C_BINDING
  USE komega_parameter, ONLY : ndim
  USE komega_vecs_c, ONLY : v3
  USE komega_math, ONLY : zcopy
  !
  IMPLICIT NONE
  !
  COMPLEX(C_DOUBLE_COMPLEX),INTENT(OUT) :: r_old(ndim)
  !
  CALL zcopy(ndim,v3,1,r_old,1)
  !
END SUBROUTINE komega_COCG_getvec
!
! Return Residual Norm
!
SUBROUTINE komega_COCG_getresidual(res) BIND(C)
  !
  USE ISO_C_BINDING
  USE komega_parameter, ONLY : nz, resnorm
  USE komega_vals_c, ONLY : pi
  !
  IMPLICIT NONE
  !
  REAL(C_DOUBLE),INTENT(OUT) :: res(nz)
  !
  res(1:nz) = resnorm / ABS(pi(1:nz))
  !
END SUBROUTINE komega_COCG_getresidual
!
! Deallocate private arrays
!
SUBROUTINE komega_COCG_finalize() BIND(C)
  !
  USE komega_parameter, ONLY : itermax, lz_conv
  USE komega_vals_c, ONLY : alpha_save, beta_save, &
  &                                 pi, pi_old, pi_save, z
  USE komega_vecs_c, ONLY : p, r_l_save, v3
  !
  IMPLICIT NONE
  !
  DEALLOCATE(z, v3, pi, pi_old, p, lz_conv)
  !
  IF(itermax > 0) THEN
     DEALLOCATE(alpha_save, beta_save, r_l_save, pi_save)
  END IF
  !
END SUBROUTINE komega_COCG_finalize
!
END MODULE komega_cocg
