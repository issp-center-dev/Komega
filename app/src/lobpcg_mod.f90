!
! This Module contains subroutines for LOBPCG calculation
!
MODULE lobpcg_mod
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC lobpcg_driver
  !
CONTAINS
!
! Driver routine for the LOBPCG
!
SUBROUTINE lobpcg_driver()
  !
  USE shiftk_vals, ONLY : ndim, rhs, e_min, e_max
  !
  IMPLICIT NONE
  !
  INTEGER :: idim, fo = 20
  COMPLEX(8),allocatable :: x(:,:), hx(:,:) ! x(:,1) = w, x(:,2) = x, x(:,3) = p
  REAL(8),allocatable :: x_r(:), x_i(:)
  !
  WRITE(*,*) 
  WRITE(*,*) "##########  Calculation of the Ground State ##########"
  WRITE(*,*)
  !
  ALLOCATE(x(ndim,3), hx(ndim,3), x_r(ndim), x_i(ndim))
  WRITE(*,*)
  WRITE(*,*) "  Compute Maximum energy"
  WRITE(*,*)
  CALL lobpcg(-1,x,hx,x_r,x_i,e_max)
  WRITE(*,*) "    E_max = ", e_max
  WRITE(*,*)
  WRITE(*,*) "  Compute Minimum energy"
  WRITE(*,*)
  CALL lobpcg(1,x,hx,x_r,x_i,e_min)
  WRITE(*,*) "    E_min = ", e_min
  !
  WRITE(*,*) 
  WRITE(*,*) "##########  Generate Right Hand Side Vector ##########"
  WRITE(*,*)
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
  OPEN(fo, file = "zvo_Excited.dat")
  WRITE(fo,*) ndim
  !
  DO idim = 1, ndim
     WRITE(fo,'(2e25.16)') DBLE(rhs(idim)), AIMAG(rhs(idim))
  END DO
  !
  CLOSE(fo)
  !
  DEALLOCATE(x, hx, x_r, x_i)
  !
END SUBROUTINE lobpcg_driver
!
! Core routine for the LOBPCG method
!
SUBROUTINE lobpcg(itarget,x,hx,x_r,x_i,eig)
  !
  USE shiftk_vals, ONLY : ndim, maxloops, threshold
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
  dnorm = SQRT(DBLE(DOT_PRODUCT(x(1:ndim,2), x(1:ndim,2))))
  x(1:ndim,2) = x(1:ndim,2) / dnorm
  CALL ham_prod(x(1:ndim,2), hx(1:ndim,2))
  x( 1:ndim,3) = CMPLX(0d0, 0d0, KIND(0d0))
  hx(1:ndim,3) = CMPLX(0d0, 0d0, KIND(0d0))
  eig = DBLE(DOT_PRODUCT(x(1:ndim,2), hx(1:ndim,2)))
  x(1:ndim, 1) = hx(1:ndim,2) - eig * x(1:ndim,2)
  !
  res = MAXVAL(ABS(x(1:ndim,1)))
  WRITE(*,'(a)')        "    iter      Residual       Energy"
  WRITE(*,'(i8,2e15.5)')        0,        res,             eig
  IF(res < threshold) GOTO 10
  !
  DO iter = 1, maxloops
     !
     CALL ham_prod(x(1:ndim,1), hx(1:ndim, 1))
     !
     DO ii = 1, 3
        DO jj = 1, 3
           hsub(jj,ii) = DOT_PRODUCT(x(1:ndim,jj), hx(1:ndim,ii)) 
           ovrp(jj,ii) = DOT_PRODUCT(x(1:ndim,jj),  x(1:ndim,ii))
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
        dnorm = SQRT(DBLE(DOT_PRODUCT(x(1:ndim,ii), x(1:ndim,ii))))
        x( 1:ndim, ii) =  x(1:ndim, ii) / dnorm
        hx(1:ndim, ii) = hx(1:ndim, ii) / dnorm
     END DO
     !
     x(1:ndim, 1) = hx(1:ndim, 2) - eig * x(1:ndim, 2)
     !
     res = MAXVAL(ABS(x(1:ndim,1)))
     WRITE(*,'(i8,2e15.5)') iter, res, eig
     IF(res < threshold) exit
     !
     dnorm = SQRT(DBLE(DOT_PRODUCT(x(1:ndim, 1), x(1:ndim, 1))))
     x(1:ndim, 1) = x(1:ndim, 1) / dnorm
     !
  END DO
  !
10 CONTINUE
  !
END SUBROUTINE lobpcg
!
END MODULE lobpcg_mod