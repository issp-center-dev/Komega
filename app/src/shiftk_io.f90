!
! Input/Output ROUTINES
!
MODULE shiftk_io
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Filname for Hamiltonian and RHS vector
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
! Read parameter for the On-The-Fly Hamiltonian generation
!
SUBROUTINE input_parameter_ham()
  !
  USE shiftk_vals, ONLY : ndim, almost0, lBiCG
  USE ham_vals, ONLY : Jx, Jy, Jz, Dz, nsite
  !
  IMPLICIT NONE
  !
  NAMELIST /ham/ Jx, Jy, Jz, Dz, nsite
  !
  Jx = 1d0
  Jy = 1d0
  Jz = 1d0
  Dz = 0d0
  nsite = 4
  !
  READ(*,ham,err=100)
  !
  ndim = 2**nsite
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Parameter for Hamiltonian ##########"
  WRITE(*,*)
  WRITE(*,*) "      Number of sites : ", nsite
  WRITE(*,*) "                   Jx : ", Jx
  WRITE(*,*) "                   Jy : ", Jy
  WRITE(*,*) "                   Jz : ", Jz
  WRITE(*,*) "                   Dz : ", Dz
  WRITE(*,*) "  Dim. of Hamiltonian : ", ndim
  !
  ! Hermitian(BiCG) or Real-Symmetric(COCG)
  !
  IF(ABS(Dz) > almost0) THEN
     WRITE(*,*) "  BiCG mathod is used."
     lBiCG = .TRUE.
  ELSE
     WRITE(*,*) "  COCG mathod is used."
     lBiCG = .FALSE.
  END IF
  !
  return
  !
100 write(*,*) "Stop in INPUT_PARAMETER for Hamiltonian. reading namelist HAM"
  !
  stop
  !
END SUBROUTINE input_parameter_ham
!
! Parameter for All CG calculations
!
SUBROUTINE input_parameter_cg()
  !
  USE shiftk_vals, ONLY : maxloops, threshold, ndim
  !
  IMPLICIT NONE
  !
  INTEGER :: convfactor
  NAMELIST /cg/ maxloops, convfactor
  !
  maxloops = ndim
  convfactor = 8
  !
  READ(*,cg,err=100)
  threshold = 10d0**(-convfactor)
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Parameter for CG Iteration ##########"
  WRITE(*,*)
  WRITE(*,*) "  Maximum number of loop : ", maxloops
  WRITE(*,*) "   Convergence Threshold : ", threshold
  !
  return
  !
100 write(*,*) "Stop in INPUT_PARAMETER for CG. reading namelist CG"
  !
  stop
  !
END SUBROUTINE input_parameter_cg
!
! Read Parameter for the Spectrum calculation
!
SUBROUTINE input_parameter_dyn()
  !
  USE shiftk_vals, ONLY : nomega, outrestart, calctype, z, e_max, e_min
  !
  IMPLICIT NONE
  !
  INTEGER :: iomega
  COMPLEX(8) :: omegamax, omegamin
  NAMELIST /dyn/ omegamax, omegamin, nomega, calctype, outrestart
  !
  calctype = "normal"
  nomega = 10
  omegamin = CMPLX(e_min, 0.01d0*(e_max - e_min), KIND(0d0))
  omegamax = CMPLX(e_max, 0.01d0*(e_max - e_min), KIND(0d0))
  outrestart = .FALSE.
  !
  READ(*,dyn,err=100)
  !
  WRITE(*,*)
  WRITE(*,*) "##########  Input Parameter for Spectrum  ##########"
  WRITE(*,*)
  WRITE(*,*) "           Max. of Omega : ", omegamax
  WRITE(*,*) "           Min. of Omega : ", omegamin
  WRITE(*,*) "           Num. of Omega : ", nomega
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
100 write(*,*) "Stop in INPUT_PARAMETER_DYN. reading namelist DYN"
  !
  stop
  !
END SUBROUTINE input_parameter_dyn
!
! Input Hamiltonian with the Matrix Market Format
!
SUBROUTINE input_hamiltonian()
  !
  USE shiftk_vals, ONLY : ndim, inham, lBiCG, almost0
  USE ham_vals, ONLY : ham, ham_indx, nham, ndiag
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
  ! Count the number of the Diagonal components
  ! Sort: Diagonal component should be first
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
  ! Hermitian(BiCG) or Real-Symmetric(COCG)
  !
  IF(MAXVAL(ABS(AIMAG(ham(1:nham)))) > almost0) THEN
     WRITE(*,*) "  BiCG mathod is used."
     lBiCG = .TRUE.
  ELSE
     WRITE(*,*) "  COCG mathod is used."
     lBiCG = .FALSE.
  END IF
  !
END SUBROUTINE input_hamiltonian
!
! Input Right Hand Side Vector
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
  ! alpha & beta for CG
  !
  DO iter = 1, iter_old
     !
     READ(fi,*) alpha_r, alpha_i, beta_r, beta_i
     alpha(iter) = CMPLX(alpha_r, alpha_i, KIND(0d0))
     beta(iter) = CMPLX(beta_r, beta_i, KIND(0d0))
     !
  END DO
  !
  ! Projected residual vector
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
  ! Last two residual vectors
  !
  DO idim = 1, ndim
     READ(fi,*) v2_r, v2_i, v12_r, v12_i
     v2(idim) = CMPLX(v2_r, v2_i, KIND(0d0))
     v12(idim) = CMPLX(v12_r, v12_i, KIND(0d0))
  END DO
  !
  ! Last two Shadow residual vectors (Only for BiCG)
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
! Output Resulting Spectrum
!
SUBROUTINE output_result()
  !
  USE shiftk_vals, ONLY : x_l, z, nomega
  !
  IMPLICIT NONE
  !
  INTEGER :: iz, fo = 20
  !
  OPEN(fo, file = "output/dynamicalG.dat")
  !
  DO iz = 1, nomega
     !
     write(fo,'(4e13.5)') DBLE(z(iz)), AIMAG(z(iz)), DBLE(x_l(1,iz)), AIMAG(x_l(1,iz))
     !
  END DO
  !
  CLOSE(fo)
  !
END SUBROUTINE output_result
!
! Output Resulting Spectrum, True residual
!
SUBROUTINE output_result_debug()
  !
  USE shiftk_vals, ONLY : v2, ndim, x_l, rhs, z, nomega
  USE ham_prod_mod, ONLY : ham_prod
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
     CALL ham_prod(x_l(1:ndim,iz), v2)
     v2(1:ndim) = z(iz) * x_l(1:ndim,iz) - v2(1:ndim) - rhs(1:ndim)
     !
     !write(*,form) v2(1:ndim)
     write(*,'(a,i5,a,2e13.5,a,e13.5)') "DEBUG (", iz, "), omega = ", z(iz), &
     &            ", Res. = ", MAXVAL(ABS(v2(1:ndim)))
     !
  END DO
  !
  OPEN(fo, file = "output/dynamicalG.dat")
  !
  DO iz = 1, nomega
     !
     Gii = dot_product(rhs,x_l(1:ndim,iz))
     write(fo,'(4e13.5)') DBLE(z(iz)), AIMAG(z(iz)), DBLE(Gii), AIMAG(Gii)
     !
  END DO
  !
  CLOSE(fo)
  !
END SUBROUTINE output_result_debug
!
END MODULE shiftk_io
