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
MODULE shifted_krylov_parameter
  !
  IMPLICIT NONE
  !
  REAL(8),PARAMETER :: &
  & almost0 = 1d-50
  !
  INTEGER,SAVE :: &
  & comm, &
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
