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
MODULE mathlib
  !
  IMPLICIT NONE
  !
  INTERFACE
     !
     REAL(8) FUNCTION ddot(n,dx,incx,dy,incy)
       REAL(8) dx(*),dy(*)
       INTEGER          incx,incy,n
     END FUNCTION ddot
     !
     COMPLEX(8) FUNCTION zdotc(n,zx,incx,zy,incy)
       COMPLEX(8) zx(*),zy(*)
       INTEGER        incx,incy,n
     END FUNCTION zdotc
     !
     SUBROUTINE  dscal(n,da,dx,incx)
       REAL(8) da,dx(*)
       INTEGER          incx,n
     END SUBROUTINE dscal
     !
     SUBROUTINE  zscal(n,za,zx,incx)
       COMPLEX(8) za,zx(*)
       INTEGER        incx,n
     END SUBROUTINE zscal
     !
     SUBROUTINE  dcopy(n,dx,incx,dy,incy)
       REAL(8) dx(*),dy(*)
       INTEGER          incx,incy,n
     END SUBROUTINE dcopy
     !
     SUBROUTINE  zcopy(n,zx,incx,zy,incy)
       COMPLEX(8) zx(*),zy(*)
       INTEGER        incx,incy,n
     END SUBROUTINE zcopy
     !
     SUBROUTINE daxpy(n,da,dx,incx,dy,incy)
       REAL(8) dx(*),dy(*),da
       INTEGER          incx,incy,n
     END SUBROUTINE daxpy
     !
     SUBROUTINE zaxpy(n,za,zx,incx,zy,incy)
       COMPLEX(8) zx(*),zy(*),za
       INTEGER        incx,incy,n
     END SUBROUTINE zaxpy
     !
     SUBROUTINE dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
       REAL(8)   alpha, beta
       INTEGER            incx, incy, lda, m, n
       CHARACTER(1)        trans
       REAL(8)   a( lda, * ), x( * ), y( * )
     END SUBROUTINE dgemv
     !
     SUBROUTINE zgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
       COMPLEX(8)         alpha, beta
       INTEGER            incx, incy, lda, m, n
       CHARACTER(1)        trans
       COMPLEX(8)         a( lda, * ), x( * ), y( * )
     END subroutine zgemv
     !
     SUBROUTINE dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
       CHARACTER(1)        transa, transb
       INTEGER            m, n, k, lda, ldb, ldc
       double precision   alpha, beta
       double precision   a( lda, * ), b( ldb, * ), c( ldc, * )
     END SUBROUTINE dgemm
     !
     SUBROUTINE zgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
       CHARACTER(1)        transa, transb
       INTEGER            m, n, k, lda, ldb, ldc
       COMPLEX(8)         alpha, beta
       COMPLEX(8)         a( lda, * ), b( ldb, * ), c( ldc, * )
     END SUBROUTINE zgemm
     !
  END INTERFACE
  !
end MODULE mathlib
