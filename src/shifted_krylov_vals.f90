MODULE shifted_krylov_parameter
  !
  IMPLICIT NONE
  !
  REAL(8),PARAMETER :: &
  & almost0 = 1d-13
  !
  INTEGER,SAVE :: &
  & iz_seed, &
  & ndim, &
  & nl,   &
  & nz,   &
  & itermax, &
  & iter
  !
  REAL(8),SAVE :: &
  & threshold
  !
END MODULE shifted_krylov_parameter
!
!
!
MODULE shifted_krylov_vals_r
  !
  IMPLICIT NONE
  !
  REAL(8),SAVE :: &
  & z_seed, &
  & rho, &
  & alpha, &
  & alpha_old, &
  & beta
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & z(:), &
  & pi(:), &
  & pi_old(:), &
  & pi_save(:,:), &
  & alpha_save(:), &
  & beta_save(:)
  !
END MODULE shifted_krylov_vals_r
!
!
!
MODULE shifted_krylov_vecs_r
  !
  IMPLICIT NONE
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & v3(:), &
  & p(:,:), &
  & r_l_save(:,:)
  !
END MODULE shifted_krylov_vecs_r
!
!
!
MODULE shifted_krylov_vals_c
  !
  IMPLICIT NONE
  !
  COMPLEX(8),SAVE :: &
  & z_seed, &
  & rho, &
  & alpha, &
  & alpha_old, &
  & beta
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & z(:), &
  & pi(:), &
  & pi_old(:), &
  & pi_save(:,:), &
  & alpha_save(:), &
  & beta_save(:)
  !
END MODULE shifted_krylov_vals_c
!
!
!
MODULE shifted_krylov_vecs_c
  !
  IMPLICIT NONE
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & v3(:), &
  & v5(:), &
  & p(:,:), &
  & r_l_save(:,:)
  !
END MODULE shifted_krylov_vecs_c
