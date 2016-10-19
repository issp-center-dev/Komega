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
