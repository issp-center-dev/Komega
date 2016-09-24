PROGRAM shiftk
  !
  USE shiftk_io, ONLY : input_filename, input_hamiltonian, input_rhs_vector, &
  &                     input_parameter_cg, input_parameter_dyn, input_parameter_ham
  USE shiftk_vals, ONLY : invec, inham
  USE lobpcg_mod, ONLY : lobpcg_driver
  USE ham_prod_mod, ONLY : print_ham, onthefly_init, finalize_ham
  USE dyn_mod, ONLY : dyn
  !
  IMPLICIT NONE
  !
  CALL system("mkdir -p output")
  !
  CALL input_filename()
  CALL input_parameter_cg()
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
  WRITE(*,*)
  WRITE(*,*) "#####  Done  #####"
  WRITE(*,*)
  !
END PROGRAM shiftk
