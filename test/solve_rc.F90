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
MODULE solve_rc_vals
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & rnd_seed, &
  & ndim,    & ! Size of Hilvert space
  & nz,      & ! Number of frequencies
  & nl,      & ! Number of Left vector
  & itermax, & ! Max. number of iteraction
  & iter_old   ! Number of iteraction of previous run
  !
  REAL(8),SAVE :: &
  & threshold ! Convergence Threshold
  !
  COMPLEX(8),SAVE :: &
  & z_seed ! Seed frequency
  !
  LOGICAL,SAVE :: &
  & restart !< It is restart run or not
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & z(:)         ! (nz): Frequency
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & ham(:,:), &
  & rhs(:), &
  & v12(:), v2(:), & ! (ndim): Working vector
  & r_l(:), & ! (nl) : Projeccted residual vector 
  & x(:,:) ! (nl,nz) : Projected result 
  !
  ! Variables for Restart
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & alpha(:), beta(:) ! (iter_old) 
  !
  COMPLEX(8),ALLOCATABLE,SAVE :: &
  & r_l_save(:,:) ! (nl,iter_old) Projected residual vectors
  !
END MODULE solve_rc_vals
!
! Routines
!
MODULE solve_rc_routines
  !
  IMPLICIT NONE
  !
CONTAINS
  !
SUBROUTINE input_size()
  !
  USE solve_rc_vals, ONLY : ndim, nl, nz, itermax, threshold, rnd_seed, restart
  !
  IMPLICIT NONE
  !
  NAMELIST /input/ ndim, nz, itermax, threshold, rnd_seed, nl, restart
  !
  ndim = 5
  nl = 5
  nz = 1
  itermax = 0
  threshold = 1d-8
  rnd_seed = 1
  restart = .FALSE.
  !
  READ(*,input,err=100)
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Standard Inputs  #####"
  WRITE(*,*)
  WRITE(*,*) "   Dimension of vvector : ", ndim
  WRITE(*,*) "  Number of frequencies : ", nz
  WRITE(*,*) "  Number of left vector : ", nl
  WRITE(*,*) "       Max. iteractions : ", itermax
  WRITE(*,*) "              Threshold : ", threshold
  WRITE(*,*) "        Seed for Random : ", threshold
  WRITE(*,*) "                Restart : ", restart
  !
  return
  !
100 write(*,*) "Stop in stdin. reading namelist file"
  !
  stop
  !
END SUBROUTINE input_size
!
! Input restart variables from file
!
SUBROUTINE input_restart()
  !
  USE solve_rc_vals, ONLY : iter_old, v2, v12, alpha, beta, z_seed, r_l_save, nl, ndim
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ierr
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Check Restart File  #####"
  WRITE(*,*)
  !
  open(fi, file = 'restart.dat',status="old", action = 'read',iostat = ierr)
  !
  IF(ierr == 0) THEN
     !
     WRITE(*,*) "  Restart file is found."
     !
     READ(fi,*) iter_old
     WRITE(*,*) "  iter_old : ", iter_old
     !
     ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl,iter_old))
     READ(fi,*) z_seed
     READ(fi,*) alpha(1:iter_old)
     READ(fi,*) beta(1:iter_old)
     READ(fi,*) r_l_save(1:nl, 1:iter_old)
     READ(fi,*) v2(1:ndim)
     READ(fi,*) v12(1:ndim)
     !
     close(fi)
     !
  ELSE
     !
     WRITE(*,*) "  Restart file is NOT found."
     STOP
     !
  END IF
  !
END SUBROUTINE input_restart
!
! Generate Equations
!
SUBROUTINE generate_system()
  !
  USE solve_rc_vals, ONLY : ndim, nz, ham, rhs, z, rnd_seed
  USE mathlib, ONLY : zgemm, zcopy
  !
  IMPLICIT NONE
  !
  INTEGER :: idim, iz
  REAL(8) :: ham_r(ndim,ndim), rhs_r(ndim), ham_i(ndim,ndim), rhs_i(ndim), rnd(rnd_seed)
  COMPLEX(8) :: ham0(ndim,ndim)
  CHARACTER(100) :: cndim, form
  !
  CALL RANDOM_NUMBER(rnd(1:rnd_seed))
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Generate Linear System  #####"
  WRITE(*,*)
  !
  READ(*,*) z(1:nz)
  WRITE(*,*) "  Frequency :"
  DO iz = 1, nz
     WRITE(*,*) iz, z(iz)
  END DO
  !
  !CALL RANDOM_NUMBER(ham(1:ndim, 1:ndim))
  ham_r(1:ndim, 1:ndim) = 0d0
  ham_i(1:ndim, 1:ndim) = 0d0
  CALL RANDOM_NUMBER(ham_r(1, 1))
  CALL RANDOM_NUMBER(ham_i(1, 1))
  DO idim = 2, ndim
     CALL RANDOM_NUMBER(ham_r(idim, idim))
     CALL RANDOM_NUMBER(ham_r(idim, idim-1))
     CALL RANDOM_NUMBER(ham_i(idim, idim))
     CALL RANDOM_NUMBER(ham_i(idim, idim-1))
  END DO
  ham(1:ndim,1:ndim) = CMPLX(ham_r(1:ndim, 1:ndim), 0d0, KIND(0d0))
  !
  CALL zgemm("C", "N", ndim, ndim, ndim, CMPLX(1d0, 0d0, KIND(0d0)), ham, ndim, ham, ndim, CMPLX(0d0, 0d0, KIND(0d0)), ham0, ndim)
  CALL zcopy(ndim*ndim,ham0,1,ham,1)
  !
  CALL RANDOM_NUMBER(rhs_r(1:ndim))
  CALL RANDOM_NUMBER(rhs_i(1:ndim))
  rhs(1:ndim) = CMPLX(0d0, 0d0, KIND(0d0))
  rhs(1) = CMPLX(1d0, 1d0, KIND(0d0))
  rhs(1:ndim) = CMPLX(rhs_r(1:ndim), rhs_i(1:ndim), KIND(0d0))
  !
  WRITE(cndim,*) ndim * 2
  WRITE(form,'(a,a,a)') "(", TRIM(ADJUSTL(cndim)), "e15.5)"
  !
  WRITE(*,*) 
  WRITE(*,*) "  Right Hand Side Vector :"
  WRITE(*,form) rhs(1:ndim)
  !
  !WRITE(*,*) 
  !WRITE(*,*) "  Matrix :"
  !DO idim = 1, ndim
  !   WRITE(*,form) ham(1:ndim,idim)
  !END DO
  !
END SUBROUTINE generate_system
!
! Output variables for restart
!
SUBROUTINE output_restart()
  !
  USE solve_rc_vals, ONLY : iter_old, v2, v12, alpha, beta, z_seed, r_l_save, nl, ndim
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Output Restart File  #####"
  WRITE(*,*)
  !
  OPEN(fo, file = 'restart.dat', action = 'write')
  !
  WRITE(fo,*) iter_old
  WRITE(fo,*) z_seed
  WRITE(fo,*) alpha(1:iter_old)
  WRITE(fo,*) beta(1:iter_old)
  WRITE(fo,*) r_l_save(1:nl, 1:iter_old)
  WRITE(fo,*) v2(1:ndim)
  WRITE(fo,*) v12(1:ndim)
  !
  close(fo)
  !
  WRITE(*,*) "  Restart File is written."
  !
END SUBROUTINE output_restart
!
! Check Result
!
SUBROUTINE output_result()
  !
  USE mathlib, ONLY : dgemv
  USE solve_rc_vals, ONLY : v2, ndim, nl, x, rhs, z, nz, ham
  !
  IMPLICIT NONE
  !
  INTEGER :: iz
  CHARACTER(100) :: cnl, form
  !
  WRITE(cnl,*) nl * 2
  WRITE(form,'(a,a,a)') "(", TRIM(ADJUSTL(cnl)), "e15.5)"
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Check Results  #####"
  WRITE(*,*)
  !
  WRITE(*,*) "  Resulting Vector"
  DO iz = 1, nz
     write(*,form) x(1:nl,iz)
  END DO
  !
  WRITE(*,*) "  Residual Vector"
  IF(nl /= ndim) THEN
     WRITE(*,*) "    Skip.  nl /= ndim."
     RETURN
  END IF
  !
  DO iz = 1, nz
     !
     v2(1:ndim) = z(iz) * x(1:ndim,iz) - rhs(1:ndim)
     CALL zgemv("N", ndim, ndim, CMPLX(-1d0, 0d0, KIND(0d0)), Ham, ndim, x(1:ndim,iz), 1, CMPLX(1d0, 0d0, KIND(0d0)), v2, 1)
     !
     write(*,form) v2(1:nl)
     !
  END DO
  !
END SUBROUTINE output_result
!
END MODULE solve_rc_routines
!
!
!
PROGRAM solve_rc
  !
  USE komega_COCG, ONLY : komega_COCG_init, komega_COCG_restart, komega_COCG_update, &
  &                       komega_COCG_getcoef, komega_COCG_getvec, komega_COCG_finalize
  USE solve_rc_routines, ONLY : input_size, input_restart, generate_system, &
  &                             output_restart, output_result
  USE solve_rc_vals, ONLY : alpha, beta, ndim, nz, nl, itermax, iter_old, ham, restart, &
  &                         rhs, v12, v2, r_l, r_l_save, threshold, x, z, z_seed
  USE mathlib, ONLY : zgemv
  !
  IMPLICIT NONE
  !
  ! Variables for Restart
  !
  INTEGER :: &
  & iter, jter, & ! Counter for Iteration
  & status(3)
  !
  COMPLEX(8),allocatable :: test_r(:,:) 
  !
  ! Input Size of vectors
  !
  CALL input_size()
  !
  ALLOCATE(v12(ndim), v2(ndim), r_l(nl), x(nl,nz), z(nz), ham(ndim,ndim), rhs(ndim))
  ALLOCATE(test_r(ndim,itermax))
  !
  CALL generate_system()
  !
  ! Check: Whether the restart file is exist.
  !
  IF (restart) CALL input_restart()
  !
  WRITE(*,*)
  WRITE(*,*) "#####  CG Initialization  #####"
  WRITE(*,*)
  !
  IF(restart) THEN
    !
    ! When restarting, counter
    !
    CALL komega_COCG_restart(ndim, nl, nz, x, z, iter_old + itermax, threshold, &
    &                 status, iter_old, v2, v12, alpha, beta, z_seed, r_l_save)
    !
    ! These vectors were saved in COCG routine
    !
    DEALLOCATE(alpha, beta, r_l_save)
    !
    IF(status(1) /= 0) GOTO 10
    !
  ELSE
     !
     ! Generate Right Hand Side Vector
     !
     v2(1:ndim) = rhs(1:ndim)
     !
     CALL komega_COCG_init(ndim, nl, nz, x, z, itermax, threshold)
     !
  END IF
  !
  ! COCG Loop
  !
  WRITE(*,*)
  WRITE(*,*) "#####  CG Iteration  #####"
  WRITE(*,*)
  !
  DO iter = 1, itermax
     !
     ! Projection of Residual vector into the space
     ! spaned by left vectors
     !
test_r(1:ndim,iter) = v2(1:ndim)
     r_l(1:nl) = v2(1:nl)
     !
     ! Matrix-vector product
     !
     CALL zgemv("N", ndim, ndim, CMPLX(1d0, 0d0, KIND(0d0)), Ham, ndim, v2, 1, CMPLX(0d0, 0d0, KIND(0d0)), v12, 1)
     !
     ! Update result x with COCG
     !
     CALL komega_COCG_update(v12, v2, x, r_l, status)
     !
     WRITE(*,'(a,i8,3i5,e15.5)') "DEBUG : ", iter, status, DBLE(v12(1))
     IF(status(1) < 0) EXIT
     !
  END DO
  !
  IF(status(2) == 0) THEN
     WRITE(*,*) "  Converged in iteration ", ABS(status(1))
  ELSE IF(status(2) == 1) THEN
     WRITE(*,*) "  Not Converged in iteration ", ABS(status(1))
  ELSE IF(status(2) == 2) THEN
     WRITE(*,*) "  Alpha becomes infinity", ABS(status(1))
  ELSE IF(status(2) == 3) THEN
     WRITE(*,*) "  Pi_seed becomes zero", ABS(status(1))
  ELSE IF(status(2) == 4) THEN
     WRITE(*,*) "  Residual & Shadow residual are orthogonal", ABS(status(1))
  END IF
  iter_old = abs(status(1))
  !
  DO iter = 1, iter_old
     DO jter = 1, iter_old
        write(*,'(e15.5)',advance="no") &
        & abs(dot_product(conjg(test_r(1:ndim,jter)), test_r(1:ndim,iter)) )
     END DO
     write(*,*)
  END DO
  !
  ! Get these vectors for restart in the Next run
  !
  ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl,iter_old))
  !
  CALL komega_COCG_getcoef(alpha, beta, z_seed, r_l_save)
  CALL komega_COCG_getvec(v12)
  !
  CALL output_restart()
  !
  DEALLOCATE(alpha, beta, r_l_save)
  !     
10 CONTINUE
  !
  ! Deallocate all intrinsic vectors
  !
  CALL komega_COCG_finalize()
  !
  ! Output to a file
  !
  CALL output_result()
  !
  DEALLOCATE(v12, v2, r_l, x, z)
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Done  #####"
  WRITE(*,*)
  !
END PROGRAM solve_rc
