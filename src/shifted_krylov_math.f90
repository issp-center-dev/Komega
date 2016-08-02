MODULE shifted_krylov_math
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
     DOUBLE COMPLEX FUNCTION zdotu(n,zx,incx,zy,incy)
       DOUBLE COMPLEX zx(*),zy(*)
       INTEGER        incx,incy,n
     END FUNCTION zdotu
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
  END INTERFACE
  !
end MODULE shifted_krylov_math
