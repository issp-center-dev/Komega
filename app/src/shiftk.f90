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
PROGRAM shiftk
  !
  USE shiftk_io, ONLY : shiftk_init, input_filename, input_hamiltonian, input_rhs_vector, &
  &                     input_parameter_cg, input_parameter_dyn, input_parameter_ham
  USE shiftk_vals, ONLY : invec, inham,stdout
  USE lobpcg_mod, ONLY : lobpcg_driver
  USE ham_prod_mod, ONLY : print_ham, onthefly_init, finalize_ham
  USE dyn_mod, ONLY : dyn
  !
  IMPLICIT NONE
  !
#if defined(MPI)
  INTEGER :: ierr
#endif
  !
  CALL shiftk_init()
  !
  CALL input_filename()
  !
  ! Initialize hamiltonian
  !
  IF(inham == "") THEN
     CALL input_parameter_ham()
     CALL onthefly_init()
#if defined(DEBUG)
     CALL print_ham()
#endif
  ELSE
     CALL input_hamiltonian()
  END IF
  !
  CALL input_parameter_cg()
  !
  ! Initialize Right Hand Side Vector
  !
  IF(invec == "") THEN
     CALL lobpcg_driver()
  ELSE
     CALL input_rhs_vector()
  END IF
  !
  ! Calculation of the Dynamical Green's function
  !
  CALL input_parameter_dyn()
  !
  CALL dyn()
  !
  CALL finalize_ham()
  !
#if defined(MPI)
  call MPI_FINALIZE(ierr)
#endif
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "#####  Done  #####"
  WRITE(stdout,*)
  !
END PROGRAM shiftk
