!
!    Copyright 2016 Mitsuaki Kawamura
!
!    This file is part of ISSP Math Library.
!
!    ISSP Math Library is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ISSP Math Library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ISSP Math Library.  If not, see <http://www.gnu.org/licenses/>.
!
!
! Routines for real-valiable CG
!
#if defined(MPI)
MODULE pshifted_bicg
#else
MODULE shifted_bicg
#endif
  !
  PRIVATE
  !
#if defined(MPI)
  PUBLIC pBiCG_init, pBiCG_restart, pBiCG_update, pBiCG_getcoef, pBiCG_getvec, pBiCG_finalize
#else
  PUBLIC BiCG_init, BiCG_restart, BiCG_update, BiCG_getcoef, BiCG_getvec, BiCG_finalize
#endif
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
  USE shifted_krylov_parameter, ONLY : iter, itermax, ndim, nz, nl, iz_seed, almost0
  USE shifted_krylov_vals_c, ONLY : alpha, alpha_save, beta_save, pi, pi_old, &
  &                               pi_save, rho, z, z_seed
  USE shifted_krylov_vecs_c, ONLY : v3, v5, r_l_save
  USE shifted_krylov_math, ONLY : dscal, zscal
  !
  IMPLICIT NONE
  !
  COMPLEX(8),INTENT(INOUT) :: v2(ndim), v4(ndim)
  INTEGER,INTENT(INOUT) :: status(3)
  !
  INTEGER :: jter
  COMPLEX(8) :: scale
  !
  status(3) = MINLOC(ABS(pi(1:nz)), 1)
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
#if defined(MPI)
SUBROUTINE pBiCG_init(ndim0, nl0, nz0, x, z0, itermax0, threshold0, comm0)
#else
SUBROUTINE BiCG_init(ndim0, nl0, nz0, x, z0, itermax0, threshold0)
#endif
  !
  USE shifted_krylov_parameter, ONLY : iter, itermax, ndim, nl, nz, &
  &                                    threshold, iz_seed
#if defined(MPI)
  USE shifted_krylov_parameter, ONLY : comm
#endif
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
#if defined(MPI)
  INTEGER,INTENT(IN) :: comm0
#endif
  !
  ndim = ndim0
  nl = nl0
  nz = nz0
  itermax = itermax0
  threshold = threshold0
#if defined(MPI)
  comm = comm0
#endif
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
  iz_seed = 1
  z_seed = z(iz_seed)
  iter = 0
  !
  IF(itermax > 0) THEN
     ALLOCATE(alpha_save(itermax), beta_save(itermax), &
     &        r_l_save(nl,itermax), pi_save(nz,-1:itermax))
     pi_save(1:nz,-1:0) = CMPLX(1d0, 0d0, KIND(0d0))
  END IF
  !
#if defined(MPI)
END SUBROUTINE pBiCG_init
#else
END SUBROUTINE BiCG_init
#endif
!
! Restart by input
!
#if defined(MPI)
SUBROUTINE pBiCG_restart(ndim0, nl0, nz0, x, z0, itermax0, threshold0, comm0, status, &
&                       iter_old, v2, v12, v4, v14, alpha_save0, beta_save0, z_seed0, r_l_save0)
#else
SUBROUTINE BiCG_restart(ndim0, nl0, nz0, x, z0, itermax0, threshold0, status, &
&                       iter_old, v2, v12, v4, v14, alpha_save0, beta_save0, z_seed0, r_l_save0)
#endif
  !
  USE shifted_krylov_parameter, ONLY : iter, itermax, ndim, nl, threshold, iz_seed
  USE shifted_krylov_vals_c, ONLY : alpha, alpha_old, alpha_save, beta, beta_save, rho, z_seed
  USE shifted_krylov_vecs_c, ONLY : r_l_save, v3, v5
  USE shifted_krylov_math, ONLY : zcopy, zdotcMPI, zabsmax
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ndim0, nl0, nz0, itermax0
  REAL(8),INTENT(IN) :: threshold0
  COMPLEX(8),INTENT(IN) :: z0(nz0)
  COMPLEX(8),INTENT(OUT) :: x(nl0,nz0)
  INTEGER,INTENT(OUT) :: status(3)
#if defined(MPI)
  INTEGER,INTENT(IN) :: comm0
#endif
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
#if defined(MPI)
  CALL pBiCG_init(ndim0, nl0, nz0, x, z0, itermax0, threshold0, comm0)
#else
  CALL BiCG_init(ndim0, nl0, nz0, x, z0, itermax0, threshold0)
#endif
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
  rho = zdotcMPI(ndim,v5,v3)
  !
  ! Seed Switching
  !
  CALL BiCG_seed_switch(v2,v4,status)
  !
  ! Convergence check
  !
  v12(1) = CMPLX(zabsmax(v2, ndim), 0d0, KIND(0d0))
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
#if defined(MPI)
END SUBROUTINE pBiCG_restart
#else
END SUBROUTINE BiCG_restart
#endif
!
! Update x, p, r
!
#if defined(MPI)
SUBROUTINE pBiCG_update(v12, v2, v14, v4, x, r_l, status)
#else
SUBROUTINE BiCG_update(v12, v2, v14, v4, x, r_l, status)
#endif
  !
  USE shifted_krylov_parameter, ONLY : iter, itermax, ndim, nl, nz, threshold, almost0
  USE shifted_krylov_vals_c, ONLY : alpha, alpha_old, alpha_save, &
  &                               beta, beta_save, rho, z_seed
  USE shifted_krylov_vecs_c, ONLY : r_l_save, v3, v5
  USE shifted_krylov_math, ONLY : zdotcMPI, zcopy, zabsmax
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
  status(1:3) = 0
  !
  rho_old = rho
  rho = zdotcMPI(ndim,v4,v2)
  IF(iter == 1) THEN
     beta = CMPLX(0d0, 0d0, KIND(0d0))
  ELSE
     beta = rho / rho_old
  END IF
  v12(1:ndim) = z_seed * v2(1:ndim) - v12(1:ndim)
  v14(1:ndim) = CONJG(z_seed) * v4(1:ndim) - v14(1:ndim)
  alpha_old = alpha
  alpha_denom = zdotcMPI(ndim,v4,v12) - beta * rho / alpha
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
  CALL BiCG_seed_switch(v2,v4,status)
  !
  ! Convergence check
  !
  v12(1) = CMPLX(zabsmax(v2, ndim), 0d0, KIND(0d0))
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
#if defined(MPI)
END SUBROUTINE pBiCG_update
#else
END SUBROUTINE BiCG_update
#endif
!
! Return saved alpha, beta, r_l
!
#if defined(MPI)
SUBROUTINE pBiCG_getcoef(alpha_save0, beta_save0, z_seed0, r_l_save0)
#else
SUBROUTINE BiCG_getcoef(alpha_save0, beta_save0, z_seed0, r_l_save0)
#endif
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
#if defined(MPI)
END SUBROUTINE pBiCG_getcoef
#else
END SUBROUTINE BiCG_getcoef
#endif
!
! Return r_old
!
#if defined(MPI)
SUBROUTINE pBiCG_getvec(r_old, r_tilde_old)
#else
SUBROUTINE BiCG_getvec(r_old, r_tilde_old)
#endif
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
#if defined(MPI)
END SUBROUTINE pBiCG_getvec
#else
END SUBROUTINE BiCG_getvec
#endif
!
! Deallocate private arrays
!
#if defined(MPI)
SUBROUTINE pBiCG_finalize()
#else
SUBROUTINE BiCG_finalize()
#endif
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
#if defined(MPI)
END SUBROUTINE pBiCG_finalize
#else
END SUBROUTINE BiCG_finalize
#endif
!
#if defined(MPI)
END MODULE pshifted_bicg
#else
END MODULE shifted_bicg
#endif
