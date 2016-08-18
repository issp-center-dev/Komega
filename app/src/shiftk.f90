PROGRAM shiftk
  !
  USE shiftk_io, ONLY : input_filename
  USE dyn_mod, ONLY : dyn
  !
  IMPLICIT NONE
  !
  CALL input_filename()
  !
  CALL dyn()
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Done  #####"
  WRITE(*,*)
  !
END PROGRAM shiftk
