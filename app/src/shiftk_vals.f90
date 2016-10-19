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
MODULE shiftk_vals
  !
  IMPLICIT NONE
  !
  REAL(8),PARAMETER :: &
  & almost0 = 1d-15
  !
  INTEGER,SAVE :: &
  & inpunit,  &
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
