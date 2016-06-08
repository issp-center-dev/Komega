MODULE syfted_krylov
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & ndim, &
  & nl,   &
  & nz,   &
  & itermax
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
  & piold(:), &
  & p(:), &
  & r(:), &
  & rold(:)
  !
CONTAINS
!
SUBROUTINE CG_R_restart(ndim, nl, nz, v2, x, z, itermax, &
    &                 threshold, r, r_old, &
    &                 alpha, beta, zseed, rsmall_save)
END SUBROUTINE CG_R_restart
!
SUBROUTINE CG_R_init(ndim0, nl0, nz0, x, z0, itermax0, threshold0)
  !
  USE mathlib, ONLY : dcopy
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: ndim, nl, nz
  REAL(8),INTENT(IN) :: x(nl), z0(nz)
  INTEGER,INTENT(IN) :: itermax0, threshold0
  !
  ndim = ndim0
  nl = nl0
  nz = nz0
  itermax = itermax0
  threshold = threshold0
  !
  ALLOCATE(z(nz), v3(ndim), pi(nz), piold(nz), p(nl))
  CALL dcopy(nz,z0,1,z,1)
  v3(1:ndim) = 0d0
  p(1:nl) = 0d0
  x(1:nl) = 0d0
  pi(1:nz) = 1d0
  piold(1:nz) = 1d0
  rho = 1d0
  alpha = 1d0
  beta = 0d0
  zseed = 0d0
  !
END SUBROUTINE CG_R_init
!
SUBROUTINE CG_R_update(v12, v2, x, rsmall, lconverged)
  !
  IMPLICIT NONE
  !
  REAL(8) :: rhoold
  !
  rhoold = rho
  rho = ddot(n,dx,incx,dy,incy)

END SUBROUTINE CG_R_update
!
SUBROUTINE CG_R_getcoef(alpha, beta, zseed, rsmall_save)
END SUBROUTINE CG_R_getcoef
!
SUBROUTINE CG_R_getvec(r, r_old)
END SUBROUTINE CG_R_getvec
!
SUBROUTINE CG_R_finalize()
  DEALLOCATE(z, v3, pi, piold, p)
END SUBROUTINE CG_R_finalize

END MODULE syfted_krylov
