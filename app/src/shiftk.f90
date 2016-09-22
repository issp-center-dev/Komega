PROGRAM shiftk
  !
  USE shiftk_io, ONLY : input_filename, input_hamiltonian, input_rhs_vector, &
  &                     input_parameter_cg, input_parameter_dyn, input_parameter_ham
  USE shiftk_vals, ONLY : invec, inham
  USE lobpcg_mod, ONLY : lobpcg_driver
  USE dyn_mod, ONLY : dyn
  !
  IMPLICIT NONE
  !
  CALL system("mkdir -p output")
  !
  CALL input_filename()
  CALL input_parameter_cg()
  !
  IF(inham == "") THEN
     CALL input_parameter_ham()
  ELSE
     CALL input_hamiltonian()
  END IF
  !
  IF(invec == "") THEN
     CALL lobpcg_driver()
  ELSE
     CALL input_rhs_vector()
  END IF
  stop
  !
  CALL input_parameter_dyn()
  !
  CALL dyn()
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Done  #####"
  WRITE(*,*)
  !
END PROGRAM shiftk
