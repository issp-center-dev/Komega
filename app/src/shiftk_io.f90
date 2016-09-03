!
! Input/Output ROUTINES
!
MODULE shiftk_io
  !
  IMPLICIT NONE
  !
CONTAINS
  !
SUBROUTINE input_filename()
  !
  USE shiftk_vals, ONLY : inham, invec
  !
  IMPLICIT NONE
  !
  NAMELIST /filename/ inham, invec
  !
  inham = "zvo_Ham.dat"
  invec = "zvo_Excited.dat"
  !
  READ(*,filename,err=100)
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input FileName  ##########"
  WRITE(*,*)
  WRITE(*,*) "  Hamiltonian file : ", TRIM(inham)
  WRITE(*,*) "  Excited state file : ", TRIM(invec)
  !
  return
  !
100 write(*,*) "Stop in INPUT_FILENAME. reading namelist FILENAME"
  !
  stop
  !
END SUBROUTINE input_filename
!
!
!
SUBROUTINE input_parameter()
  !
  USE shiftk_vals, ONLY : nomega, maxloops, outrestart, threshold, ndim, calctype, z
  !
  IMPLICIT NONE
  !
  INTEGER :: convfactor, iomega
  COMPLEX(8) :: omegamax, omegamin
  NAMELIST /parameter/ omegamax, omegamin, nomega, maxloops, calctype, convfactor, outrestart
  !
  maxloops = ndim
  calctype = "normal"
  convfactor = 8
  nomega = 10
  omegamin = CMPLX(0d0, 1d0, KIND(0d0))
  omegamax = CMPLX(1d0, 1d0, KIND(0d0))
  outrestart = .FALSE.
  !
  READ(*,parameter,err=100)
  threshold = 10d0**(-convfactor)
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Parameter  ##########"
  WRITE(*,*)
  WRITE(*,*) "           Max. of Omega : ", omegamax
  WRITE(*,*) "           Min. of Omega : ", omegamin
  WRITE(*,*) "           Num. of Omega : ", nomega
  WRITE(*,*) "  Maximum number of loop : ", maxloops
  WRITE(*,*) "   Convergence Threshold : ", threshold
  WRITE(*,'(a,a)') "         Calculation type : ", calctype
  !
  ALLOCATE(z(nomega))
  z(1) = omegamin
  DO iomega = 2, nomega
     z(iomega) = omegamin + (omegamax - omegamin) * DBLE(iomega - 1) / DBLE(nomega - 1)
  END DO
  !
  return
  !
100 write(*,*) "Stop in INPUT_PARAMETER. reading namelist PARAMETER"
  !
  stop
  !
END SUBROUTINE input_parameter
!
! Input Hamiltonian with the Matrix Market Format
!
SUBROUTINE input_hamiltonian()
  !
  USE shiftk_vals, ONLY : ndim, ham, ham_indx, nham, ndiag, inham, lBiCG
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ndim2, iham, ham_indx0(2)
  REAL(8) :: ham_r, ham_i
  COMPLEX(8) :: ham0
  CHARACTER(100) :: ctmp
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Hamiltonian  ##########"
  WRITE(*,*)
  !
  OPEN(fi, file = TRIM(inham))
  READ(fi, *) ctmp
  READ(fi,*) ndim, ndim2, nham
  WRITE(*,*) "          Dim. of Hamiltonian : ", ndim, ndim2
  WRITE(*,*) "  Num. of Non-Zero Components : ", nham
  !
  IF(ndim2 /= ndim) THEN
     WRITE(*,*) "ERROR ! Hamiltonian is not square."
     STOP
  END IF
  !
  ALLOCATE(ham(nham), ham_indx(2,nham))
  !
  DO iham = 1, nham
     READ(fi,*) ham_indx(1,iham), ham_indx(2,iham), ham_r, ham_i
     ham(iham) = CMPLX(ham_r, ham_i, KIND(0d0))
  END DO
  !
  ndiag = 0
  DO iham = 1, nham
     !
     IF(ham_indx(1,iham) == ham_indx(2,iham)) THEN
        !
        ndiag = ndiag + 1
        !
        ham_indx0(1:2) = ham_indx(1:2,iham)
        ham_indx(1:2,iham) = ham_indx(1:2,ndiag)
        ham_indx(1:2,ndiag) = ham_indx0(1:2)
        !
        ham0 = ham(iham)
        ham(iham) = ham(ndiag)
        ham(ndiag) = ham0
        !
     END IF
     !
  END DO
  !
  CLOSE(fi)
  !
  IF(MAXVAL(ABS(AIMAG(ham(1:nham)))) > 1d-8) THEN
     WRITE(*,*) "  BiCG mathod is used."
     lBiCG = .TRUE.
  ELSE
     WRITE(*,*) "  COCG mathod is used."
     lBiCG = .FALSE.
  END IF
  !
END SUBROUTINE input_hamiltonian
!
! Input Excited State
!
SUBROUTINE input_rhs_vector()
  !
  USE shiftk_vals, ONLY : ndim, rhs, invec
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ndim2, idim
  REAL(8) :: rhs_r, rhs_i
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Right Hand Side Vector ##########"
  WRITE(*,*)
  OPEN(fi, file = TRIM(invec))
  READ(fi,*) ndim2
  WRITE(*,*) "  Dim. of RHS vector : ", ndim2
  !
  IF(ndim2 /= ndim) THEN
     WRITE(*,*) "ERROR ! Dimension is Incorrect."
     STOP
  END IF
  !
  ALLOCATE(rhs(ndim))
  !
  DO idim = 1, ndim
     READ(fi,*) rhs_r, rhs_i
     rhs(idim) = CMPLX(rhs_r, rhs_i, KIND(0d0))
  END DO
  !
  CLOSE(fi)
  !
END SUBROUTINE input_rhs_vector
!
! Input restart variables from file
!
SUBROUTINE input_restart_parameter()
  !
  USE shiftk_vals, ONLY : iter_old, alpha, beta, z_seed, nl, r_l_save
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, iter, il
  REAL(8) :: z_seed_r, z_seed_i, alpha_r, alpha_i, beta_r, beta_i
  REAL(8) :: r_l_save_r, r_l_save_i
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Restart Parameter  ##########"
  WRITE(*,*)
  !
  OPEN(fi, file = 'output/TriDiagComp.dat')
  !
  READ(fi,*) iter_old
  WRITE(*,*) "  Num. of Iteration (Previous Run) : ", iter_old
  !
  ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl, iter_old))
  !
  READ(fi,*) z_seed_r, z_seed_i
  WRITE(*,*) "  Previous Omega_Seed : ", z_seed_r, z_seed_i
  z_seed = CMPLX(z_seed_r, z_seed_i, KIND(0d0))
  !
  DO iter = 1, iter_old
     !
     READ(fi,*) alpha_r, alpha_i, beta_r, beta_i
     alpha(iter) = CMPLX(alpha_r, alpha_i, KIND(0d0))
     beta(iter) = CMPLX(beta_r, beta_i, KIND(0d0))
     !
  END DO
  !
  DO iter = 1, iter_old
     !
     DO il = 1, nl
        READ(fi,*) r_l_save_r, r_l_save_i
        r_l_save(il,iter) = CMPLX(r_l_save_r, r_l_save_i, KIND(0d0))
     END DO
     !
  END DO
  !
  CLOSE(fi)
  !
END SUBROUTINE input_restart_parameter
!
! Input restart variables from file
!
SUBROUTINE input_restart_vector()
  !
  USE shiftk_vals, ONLY : ndim, v2, v12, v4, v14, lBiCG
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ndim2, idim
  REAL(8) :: v2_r, v2_i, v12_r, v12_i
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Restart Vector  ##########"
  WRITE(*,*)
  !
  OPEN(fi, file = 'output/ResVec.dat')
  !
  READ(fi,*) ndim2
  WRITE(*,*) "  Dim. of Residual vector : ", ndim2
  !
  IF(ndim2 /= ndim) THEN
     WRITE(*,*) "ERROR ! Dimension is Incorrect."
     STOP
  END IF
  !
  DO idim = 1, ndim
     READ(fi,*) v2_r, v2_i, v12_r, v12_i
     v2(idim) = CMPLX(v2_r, v2_i, KIND(0d0))
     v12(idim) = CMPLX(v12_r, v12_i, KIND(0d0))
  END DO
  !
  IF(lBiCG) THEN
     READ(fi,*) v2_r, v2_i, v12_r, v12_i
     v4(idim) = CMPLX(v2_r, v2_i, KIND(0d0))
     v14(idim) = CMPLX(v12_r, v12_i, KIND(0d0))
  END IF
  !
  CLOSE(fi)
  !
END SUBROUTINE input_restart_vector
!
! Input restart variables from file
!
SUBROUTINE output_restart_parameter()
  !
  USE shiftk_vals, ONLY : iter_old, alpha, beta, z_seed, r_l_save, nl
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, iter, il
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Output Restart Parameter  ##########"
  WRITE(*,*)
  !
  OPEN(fo, file = 'output/TriDiagComp.dat')
  !
  WRITE(fo,*) iter_old
  WRITE(*,*) "  Num. of Iteration (Current Run) : ", iter_old
  !
  WRITE(*,'(a,2e13.5)') "   Current Omega_Seed : ", z_seed
  WRITE(fo,'(2e25.16)') z_seed
  !
  DO iter = 1, iter_old
     !
     WRITE(fo,'(4e25.16)') alpha(iter), beta(iter)
     !
  END DO
  !
  DO iter = 1, iter_old
     !
     DO il = 1, nl
        WRITE(fo,'(2e25.16)') r_l_save(il,iter)
     END DO
     !
  END DO
  !
  CLOSE(fo)
  !
END SUBROUTINE output_restart_parameter
!
! Input restart variables from file
!
SUBROUTINE output_restart_vector()
  !
  USE shiftk_vals, ONLY : ndim, v2, v12, v4, v14, lBiCG
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, idim
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Output Restart Vector  ##########"
  WRITE(*,*)
  !
  OPEN(fo, file = 'output/ResVec.dat')
  !
  WRITE(fo,*) ndim
  WRITE(*,*) "  Dim. of Residual vector : ", ndim
  !
  DO idim = 1, ndim
     WRITE(fo,'(4e25.16)') v2(idim), v12(idim)
  END DO
  !
  IF(lBiCG) THEN
     DO idim = 1, ndim
        WRITE(fo,'(4e25.16)') v4(idim), v14(idim)
     END DO
  END IF
  !
  CLOSE(fo)
  !
END SUBROUTINE output_restart_vector
!
! Check Result
!
SUBROUTINE output_result()
  !
  USE shiftk_vals, ONLY : x, z, nomega
  !
  IMPLICIT NONE
  !
  INTEGER :: iz, fo = 20
  !
  OPEN(fo, file = "output/dynamicalG.dat")
  !
  DO iz = 1, nomega
     !
     write(fo,'(4e13.5)') DBLE(z(iz)), AIMAG(z(iz)), DBLE(x(1,iz)), AIMAG(x(1,iz))
     !
  END DO
  !
  CLOSE(fo)
  !
END SUBROUTINE output_result
!
! Check Result
!
SUBROUTINE output_result_debug()
  !
  USE shiftk_vals, ONLY : v2, ndim, x, rhs, z, nomega
  USE ham_prod, ONLY : ham_prod_compress
  !
  IMPLICIT NONE
  !
  INTEGER :: iz, fo = 20
  CHARACTER(100) :: cndim, form
  COMPLEX(8) :: Gii
  !
  WRITE(cndim,*) ndim * 2
  WRITE(form,'(a,a,a)') "(", TRIM(ADJUSTL(cndim)), "e13.5)"
  !
  WRITE(*,*)
  WRITE(*,*) "#####  Check Results  #####"
  WRITE(*,*)
  !
  WRITE(*,*) "  Residual Vector"
  !
  DO iz = 1, nomega
     !
     CALL ham_prod_compress(x(1:ndim,iz), v2)
     v2(1:ndim) = z(iz) * x(1:ndim,iz) - v2(1:ndim) - rhs(1:ndim)
     !
     !write(*,form) v2(1:ndim)
     write(*,'(a,i5,a,2e13.5,a,e13.5)') "DEBUG (", iz, "), omega = ", z(iz), &
     &            ", Res. = ", dble(dot_product(v2,v2))
     !
  END DO
  !
  OPEN(fo, file = "output/dynamicalG.dat")
  !
  DO iz = 1, nomega
     !
     Gii = dot_product(rhs,x(1:ndim,iz))
     write(fo,'(4e13.5)') DBLE(z(iz)), AIMAG(z(iz)), DBLE(Gii), AIMAG(Gii)
     !
  END DO
  !
  CLOSE(fo)
  !
END SUBROUTINE output_result_debug
!
END MODULE shiftk_io
