MODULE shiftk_vals
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & ndiag,   & ! Diagonal components
  & nham,    & ! Non-zero elements of compressed Hamiltonian
  & ndim,    & ! Size of Hilvert space
  & nl,      & ! Dimention of x
  & nomega,      & ! Number of frequencies
  & maxloops, & ! Max. number of iteraction
  & iter_old   ! Number of iteraction of previous run
  !
  REAL(8),SAVE :: &
  & threshold ! Convergence Threshold
  !
  COMPLEX(8),SAVE :: &
  & z_seed ! Seed frequency
  !
  INTEGER,ALLOCATABLE,SAVE :: &
  & ham_indx(:,:) ! row & column index of Hamiltonian
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & z(:)         ! (nomega): Frequency
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & ham(:), & ! Compressed Hamiltonian
  & rhs(:), &
  & v12(:), v2(:), & ! (ndim): Working vector
  & v14(:), v4(:), & ! (ndim): Working vector
  & r_l(:), & ! (nl) : Projeccted residual vector 
  & x(:,:) ! (nl,nomega) : Projected result 
  !
  ! Variables for Restart
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & alpha(:), beta(:) ! (iter_old) 
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & r_l_save(:,:) ! (nl,iter_old) Projected residual vectors
  !
  CHARACTER(256),SAVE :: &
  & inham, &
  & invec, &
  & calctype
  !
  LOGICAL,SAVE :: &
  lBiCG, &
  outrestart
  !
END MODULE shiftk_vals
