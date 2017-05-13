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
