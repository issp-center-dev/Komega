MODULE shiftk_vals
  !
  IMPLICIT NONE
  !
  REAL(8),PARAMETER :: &
  & almost0 = 1d-15
  !
  INTEGER,SAVE :: &
  & stdout,   &
  & nproc,    &
  & myrank,  &
  & ndim,     & ! Size of Hilvert space
  & nl,       & ! Dimention of x
  & nomega,   & ! Number of frequencies
  & maxloops, & ! Max. number of iteraction
  & iter_old   ! Number of iteraction of previous run
  !
  REAL(8),SAVE :: &
  & e_min,  & ! Minimum energy
  & e_max,  & ! Maximum energy
  & threshold ! Convergence Threshold
  !
  COMPLEX(8),SAVE :: &
  & z_seed ! Seed frequency
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & z(:) ! (nomega): Frequency
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & rhs(:),        &
  & v12(:), v2(:), & ! (ndim): Working vector
  & v14(:), v4(:), & ! (ndim): Working vector
  & r_l(:),        & ! (nl) : Projeccted residual vector 
  & x_l(:,:)         ! (nl,nomega) : Projected result 
  !
  CHARACTER(256),SAVE :: &
  & inham, & ! File name for Hamiltonian
  & invec, & ! File name for RHS vector
  & calctype ! Restart type
  !
  LOGICAL,SAVE :: &
  lBiCG, &   ! BiCG is used for Complex-Hermitian case
  outrestart ! Flag for output restart parameters
  !
  ! Variables for Restart
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & alpha(:), & ! (iter_old) CG alpha
  & beta(:)     ! (iter_old) CG beta
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & r_l_save(:,:) ! (nl,iter_old) Projected residual vectors
  !
END MODULE shiftk_vals
