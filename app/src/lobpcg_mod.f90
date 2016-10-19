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
! This Module contains subroutines for LOBPCG calculation
!
MODULE lobpcg_mod
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC lobpcg_driver, zabsmax, zdotcMPI
  !
  INTERFACE
     DOUBLE COMPLEX FUNCTION zdotc(n,zx,incx,zy,incy)
       DOUBLE COMPLEX zx(*),zy(*)
       INTEGER        incx,incy,n
     END FUNCTION zdotc
  END INTERFACE
  !
CONTAINS
!
! Driver routine for the LOBPCG
!
SUBROUTINE lobpcg_driver()
  !
  USE shiftk_vals, ONLY : ndim, rhs, e_min, e_max, stdout
#if defined(DEBUG)
  USE shiftk_vals, ONLY : nproc
#endif  
  !
  IMPLICIT NONE
  !
  INTEGER :: idim
#if defined(DEBUG)
  INTEGER :: fo = 20
#endif
  COMPLEX(8),allocatable :: x(:,:), hx(:,:) ! x(:,1) = w, x(:,2) = x, x(:,3) = p
  REAL(8),allocatable :: x_r(:), x_i(:)
  !
  WRITE(stdout,*) 
  WRITE(stdout,*) "##########  Calculation of the Ground State ##########"
  WRITE(stdout,*)
  !
  ALLOCATE(x(ndim,3), hx(ndim,3), x_r(ndim), x_i(ndim))
  WRITE(stdout,*)
  WRITE(stdout,*) "  Compute Maximum energy"
  WRITE(stdout,*)
  CALL lobpcg(-1,x,hx,x_r,x_i,e_max)
  WRITE(stdout,*) "    E_max = ", e_max
  WRITE(stdout,*)
  WRITE(stdout,*) "  Compute Minimum energy"
  WRITE(stdout,*)
  CALL lobpcg(1,x,hx,x_r,x_i,e_min)
  WRITE(stdout,*) "    E_min = ", e_min
  !
  WRITE(stdout,*) 
  WRITE(stdout,*) "##########  Generate Right Hand Side Vector ##########"
  WRITE(stdout,*)
  !
  ALLOCATE(rhs(ndim))
  !
  DO idim = 1, ndim
     IF(MOD(idim, 2) == 1) THEN
        rhs(idim) = - 0.5d0 * x(idim,2)
     ELSE 
        rhs(idim) =   0.5d0 * x(idim,2)
     END IF
  END DO
  !
#if defined(DEBUG)
  IF(nproc == 1) THEN
     !
     OPEN(fo, file = "zvo_Excited.dat")
     WRITE(fo,*) ndim
     !
     DO idim = 1, ndim
        WRITE(fo,'(2e25.16)') DBLE(rhs(idim)), AIMAG(rhs(idim))
     END DO
     !
     CLOSE(fo)
     !
  END IF
#endif
  !
  DEALLOCATE(x, hx, x_r, x_i)
  !
END SUBROUTINE lobpcg_driver
!
! Core routine for the LOBPCG method
!
SUBROUTINE lobpcg(itarget,x,hx,x_r,x_i,eig)
  !
  USE shiftk_vals, ONLY : ndim, maxloops, threshold, stdout
  USE ham_prod_mod, ONLY : ham_prod
  !  
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: itarget
  COMPLEX(8),INTENT(OUT) :: x(ndim,3), hx(ndim,3) ! x(:,1) = w, x(:,2) = x, x(:,3) = p
  REAL(8),INTENT(OUT) :: x_r(ndim), x_i(ndim)
  REAL(8),INTENT(OUT) :: eig
  !
  INTEGER :: lwork = 5, info, iter, ii, jj, jtarget
  REAL(8) :: dnorm, rwork(7), eig3(3), res
  COMPLEX(8) :: hsub(3,3), ovrp(3,3), work(5)
  !
  ! Initial guess of re
  !
  CALL RANDOM_NUMBER(x_r(1:ndim))
  CALL RANDOM_NUMBER(x_i(1:ndim))
  x(1:ndim,2) = CMPLX(x_r(1:ndim), x_i(1:ndim), KIND(0d0))
  dnorm = SQRT(DBLE(zdotcMPI(ndim, x(1:ndim,2), x(1:ndim,2))))
  x(1:ndim,2) = x(1:ndim,2) / dnorm
  CALL ham_prod(x(1:ndim,2), hx(1:ndim,2))
  x( 1:ndim,3) = CMPLX(0d0, 0d0, KIND(0d0))
  hx(1:ndim,3) = CMPLX(0d0, 0d0, KIND(0d0))
  eig = DBLE(zdotcMPI(ndim, x(1:ndim,2), hx(1:ndim,2)))
  x(1:ndim, 1) = hx(1:ndim,2) - eig * x(1:ndim,2)
  !
  res = zabsmax(x(1:ndim,1), ndim)
  WRITE(stdout,'(a)')        "    iter      Residual       Energy"
  WRITE(stdout,'(i8,2e15.5)')        0,        res,             eig
  IF(res < threshold) GOTO 10
  !
  DO iter = 1, maxloops
     !
     CALL ham_prod(x(1:ndim,1), hx(1:ndim, 1))
     !
     DO ii = 1, 3
        DO jj = 1, 3
           hsub(jj,ii) = zdotcMPI(ndim, x(1:ndim,jj), hx(1:ndim,ii)) 
           ovrp(jj,ii) = zdotcMPI(ndim, x(1:ndim,jj),  x(1:ndim,ii))
        END DO
     END DO
     eig = DBLE(hsub(2,2))
     !
     IF(iter == 1) THEN
        CALL zhegv(1, 'V', 'U', 2, hsub, 3, ovrp, 3, eig3, work, lwork, rwork, info)
        eig3(3) = 0d0
        !
        IF(itarget == 1) THEN
           jtarget = 1
        ELSE
           jtarget = 2
        END IF
        !
     ELSE
        CALL zhegv(1, 'V', 'U', 3, hsub, 3, ovrp, 3, eig3, work, lwork, rwork, info)
        !
        IF(itarget == 1) THEN
           jtarget = 1
        ELSE
           jtarget = 3
        END IF
        !
    END IF
     !
     eig = 0.5d0 * (eig + eig3(jtarget))
     !
     x( 1:ndim,2) = MATMUL( x(1:ndim,1:3), hsub(1:3,jtarget))
     hx(1:ndim,2) = MATMUL(hx(1:ndim,1:3), hsub(1:3,jtarget))
     x( 1:ndim,3) = hsub(1,jtarget) *  x(1:ndim,1) + hsub(3,jtarget) *  x(1:ndim,3)
     hx(1:ndim,3) = hsub(1,jtarget) * hx(1:ndim,1) + hsub(3,jtarget) * hx(1:ndim,3)
     !
     DO ii = 2, 3
        dnorm = SQRT(DBLE(zdotcMPI(ndim, x(1:ndim,ii), x(1:ndim,ii))))
        x( 1:ndim, ii) =  x(1:ndim, ii) / dnorm
        hx(1:ndim, ii) = hx(1:ndim, ii) / dnorm
     END DO
     !
     x(1:ndim, 1) = hx(1:ndim, 2) - eig * x(1:ndim, 2)
     !
     res = zabsmax(x(1:ndim,1), ndim)
     WRITE(stdout,'(i8,2e15.5)') iter, res, eig
     IF(res < threshold) exit
     !
     dnorm = SQRT(DBLE(zdotcMPI(ndim, x(1:ndim, 1), x(1:ndim, 1))))
     x(1:ndim, 1) = x(1:ndim, 1) / dnorm
     !
  END DO
  !
10 CONTINUE
  !
END SUBROUTINE lobpcg
!
! MAXVAL with MPI allreduce (for complex(8))
!
FUNCTION zabsmax(array, n) RESULT(maxarray)
  !
#if defined(MPI)
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD
#endif
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n
  COMPLEX(8),INTENT(IN) :: array(n)
  REAL(8) maxarray
  !
#if defined(MPI)
  INTEGER :: ierr
#endif
  !
  maxarray = MAXVAL(ABS(array))
  !
#if defined(MPI)
  call MPI_allREDUCE(MPI_IN_PLACE, maxarray, 1, &
  &                  MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
  !
END FUNCTION zabsmax
!
! zdotc with MPI allreduce
!
FUNCTION zdotcMPI(n,zx,zy) RESULT(prod)
  !
#if defined(MPI)
  use mpi, only : MPI_IN_PLACE, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD
#endif
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(IN) :: n
  COMPLEX(8),INTENT(IN) :: zx(n), zy(n)
  COMPLEX(8) prod
  !
#if defined(MPI)
  INTEGER :: ierr
#endif
  !
#if defined(NO_ZDOT)
  prod = DOT_PRODUCT(zx,zy)
#else
  prod = zdotc(n,zx,1,zy,1)
#endif
  !
#if defined(MPI)
  call MPI_allREDUCE(MPI_IN_PLACE, prod, 1, &
  &                  MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  !
END FUNCTION zdotcMPI
!
END MODULE lobpcg_mod
