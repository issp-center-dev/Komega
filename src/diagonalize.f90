MODULE variables
  !
  IMPLICIT NONE
  !
  INTEGER,SAVE :: &
  & ndim, &
  & ninter
  !
  REAL(8),ALLOCATABLE,SAVE :: &
  & diag(:)
  !
END MODULE variables

SUBROUTINE matrix_mpultiply(vin, vout)
  !
  USE variables, ONLY : ndim, diag
  !
  IMPLICIT NONE
  !
  REAL(8),INTENT(IN) :: vin(ndim)
  REAL(8),INTENT(OUT) :: vout(ndim)
  !
  vout(1:ndim) = diag(1:ndim) * vin(1:ndim)
  !
  ! Exchange interaction
  !
  DO intr = 1, nexc0

  END DO

	  sigma1=0; sigma2=1;

	    child_exchange_spin_GetInfo(i, X);
	    dam_pr = GC_child_exchange_spin(tmp_v0, tmp_v1, X);
          }
	  X->Large.prdct += dam_pr;
	}/* for (i = 0; i < X->Def.NExchangeCoupling; i ++) */




END SUBROUTINE matrix_mpultiply



PROGRAM diagonalize
  !
  IMPLICIT NONE
  !
END PROGRAM diagonalize
