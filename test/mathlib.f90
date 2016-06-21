MODULE mathlib
  !
  IMPLICIT NONE
  !
  INTERFACE
     !
     DOUBLE PRECISION FUNCTION ddot(n,dx,incx,dy,incy)
       DOUBLE PRECISION dx(*),dy(*)
       INTEGER          incx,incy,n
     END FUNCTION ddot
     !
     DOUBLE COMPLEX FUNCTION zdotc(n,zx,incx,zy,incy)
       DOUBLE COMPLEX zx(*),zy(*)
       INTEGER        incx,incy,n
     END FUNCTION zdotc
     !
     SUBROUTINE  dscal(n,da,dx,incx)
       DOUBLE PRECISION da,dx(*)
       INTEGER          incx,n
     END SUBROUTINE dscal
     !
     SUBROUTINE  zscal(n,za,zx,incx)
       DOUBLE COMPLEX za,zx(*)
       INTEGER        incx,n
     END SUBROUTINE zscal
     !
     SUBROUTINE  dcopy(n,dx,incx,dy,incy)
       DOUBLE PRECISION dx(*),dy(*)
       INTEGER          incx,incy,n
     END SUBROUTINE dcopy
     !
     SUBROUTINE  zcopy(n,zx,incx,zy,incy)
       DOUBLE COMPLEX zx(*),zy(*)
       INTEGER        incx,incy,n
     END SUBROUTINE zcopy
     !
     SUBROUTINE daxpy(n,da,dx,incx,dy,incy)
       DOUBLE PRECISION dx(*),dy(*),da
       INTEGER          incx,incy,n
     END SUBROUTINE daxpy
     !
     SUBROUTINE zaxpy(n,za,zx,incx,zy,incy)
       DOUBLE COMPLEX zx(*),zy(*),za
       INTEGER        incx,incy,n
     END SUBROUTINE zaxpy
     !
     SUBROUTINE dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
       DOUBLE PRECISION   alpha, beta
       INTEGER            incx, incy, lda, m, n
       CHARACTER*1        trans
       DOUBLE PRECISION   a( lda, * ), x( * ), y( * )
     END SUBROUTINE dgemv
     !
     SUBROUTINE dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
       CHARACTER*1        transa, transb
       INTEGER            m, n, k, lda, ldb, ldc
       double precision   alpha, beta
       double precision   a( lda, * ), b( ldb, * ), c( ldc, * )
     END SUBROUTINE dgemm
     !
     SUBROUTINE zgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
       CHARACTER*1        transa, transb
       INTEGER            m, n, k, lda, ldb, ldc
       COMPLEX*16         alpha, beta
       COMPLEX*16         a( lda, * ), b( ldb, * ), c( ldc, * )
     END SUBROUTINE zgemm
     !
  END INTERFACE
  !
end MODULE mathlib
