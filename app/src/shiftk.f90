PROGRAM shiftk
  !
  USE shiftk_io, ONLY : input_filename
  USE dyn_mod, ONLY : dyn
  !
  IMPLICIT NONE
  !
  CALL system("mkdir -p output")
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
