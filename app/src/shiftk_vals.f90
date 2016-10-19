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
