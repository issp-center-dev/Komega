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
  & z_seed, &
  & rho, &
  & alpha, &
  & alpha_old, &
  & beta
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & z(:), &
  & v3(:), &
  & pi(:), &
  & pi_old(:), &
  & pi_save(:,:), &
  & p(:,:), &
  & alpha_save(:), &
  & beta_save(:), &
  & r_l_save(:,:)
  !
END MODULE shifted_krylov_vals
