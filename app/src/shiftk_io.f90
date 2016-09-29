!
! Input/Output ROUTINES
!
MODULE shiftk_io
  !
  IMPLICIT NONE
  !
CONTAINS
!
! Initialization of MPI
!
SUBROUTINE shiftk_init()
  !
#if defined(MPI)
  USE mpi, only : MPI_COMM_WORLD
#endif
  USE shiftk_vals, ONLY : myrank, nproc, stdout
  !
  IMPLICIT NONE
  !
#if defined(MPI)
  INTEGER ierr
#endif
  !  
#if defined(MPI)
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, nproc, ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, myrank, ierr)
#else
  nproc = 1
  myrank = 0
#endif
  !
  IF(myrank == 0) THEN
     CALL system("mkdir -p output")
     stdout = 6
  ELSE
     stdout = 6
     OPEN(stdout, file='/dev/null', status='unknown')
  END IF   
  !
END SUBROUTINE shiftk_init
!
! Filname for Hamiltonian and RHS vector
!
SUBROUTINE input_filename()
  !
#if defined(MPI)
  USE mpi, only : MPI_COMM_WORLD, MPI_CHARACTER
#endif
  USE shiftk_vals, ONLY : inham, invec, stdout, myrank
  !
  IMPLICIT NONE
  !
#if defined(MPI)
  INTEGER ierr
#endif
  NAMELIST /filename/ inham, invec
  !
  inham = ""
  invec = ""
  !
  IF(myrank == 0) READ(*,filename,err=100)
  !
#if defined(MPI)
  call MPI_BCAST(inham, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(invec, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
#endif
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input FileName  ##########"
  WRITE(stdout,*)
  WRITE(stdout,*) "  Hamiltonian file : ", TRIM(inham)
  WRITE(stdout,*) "  Excited state file : ", TRIM(invec)
  !
  return
  !
100 write(*,*) "Stop in INPUT_FILENAME. reading namelist FILENAME"
  !
#if defined(MPI)
     CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  STOP
  !
END SUBROUTINE input_filename
!
! Read parameter for the On-The-Fly Hamiltonian generation
!
SUBROUTINE input_parameter_ham()
  !
#if defined(MPI)
  USE mpi, only : MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_PRECISION
#endif
  USE shiftk_vals, ONLY : ndim, almost0, lBiCG, stdout, myrank, nproc
  USE ham_vals, ONLY : Jx, Jy, Jz, Dz, nsite, nsitep
  !
  IMPLICIT NONE
  !
#if defined(MPI)
  INTEGER ierr
#endif
  NAMELIST /ham/ Jx, Jy, Jz, Dz, nsite
  !
  Jx = 1d0
  Jy = 1d0
  Jz = 1d0
  Dz = 0d0
  nsite = 4
  !
  IF(myrank == 0) READ(*,ham,err=100)
  !
#if defined(MPI)
  call MPI_BCAST(nsite, 1, MPI_INTEGER,       0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(Jx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(Jy, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(Jz, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(Dz, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Parameter for Hamiltonian ##########"
  WRITE(stdout,*)
  WRITE(stdout,*) "    TOTAL Number of sites : ", nsite
  !
  nsitep = NINT(LOG(DBLE(nproc)) / LOG(2d0))
  IF(2**nsitep /= nproc) THEN
     WRITE(*,*) "ERROR ! Number of processes is not 2-exponent."
#if defined(MPI)
     CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
     STOP
  END IF
  WRITE(stdout,*) "    LOCAL Number of sites : ", nsite - nsitep
  !
  WRITE(stdout,*) "                       Jx : ", Jx
  WRITE(stdout,*) "                       Jy : ", Jy
  WRITE(stdout,*) "                       Jz : ", Jz
  WRITE(stdout,*) "                       Dz : ", Dz
  ndim = 2**(nsite - nsitep)
  WRITE(stdout,*) "  Dim. of Hamiltonian : ", ndim
  !
  ! Hermitian(BiCG) or Real-Symmetric(COCG)
  !
  IF(ABS(Dz) > almost0) THEN
     WRITE(stdout,*) "  BiCG mathod is used."
     lBiCG = .TRUE.
  ELSE
     WRITE(stdout,*) "  COCG mathod is used."
     lBiCG = .FALSE.
  END IF
  !
  return
  !
100 write(*,*) "Stop in INPUT_PARAMETER for Hamiltonian. reading namelist HAM"
  !
#if defined(MPI)
  CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  STOP
  !
END SUBROUTINE input_parameter_ham
!
! Parameter for All CG calculations
!
SUBROUTINE input_parameter_cg()
  !
#if defined(MPI)
  USE mpi, only : MPI_COMM_WORLD, MPI_INTEGER
#endif
  USE shiftk_vals, ONLY : maxloops, threshold, ndim, stdout, myrank
  !
  IMPLICIT NONE
  !
  INTEGER :: convfactor
#if defined(MPI)
  INTEGER ierr
#endif
  NAMELIST /cg/ maxloops, convfactor
  !
  maxloops = ndim
  convfactor = 8
  !
  IF(myrank == 0) READ(*,cg,err=100)
  !
#if defined(MPI)
  call MPI_BCAST(maxloops,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(convfactor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
  !
  threshold = 10d0**(-convfactor)
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Parameter for CG Iteration ##########"
  WRITE(stdout,*)
  WRITE(stdout,*) "  Maximum number of loop : ", maxloops
  WRITE(stdout,*) "   Convergence Threshold : ", threshold
  !
  return
  !
100 write(*,*) "Stop in INPUT_PARAMETER for CG. reading namelist CG"
  !
#if defined(MPI)
  CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  STOP
  !
END SUBROUTINE input_parameter_cg
!
! Read Parameter for the Spectrum calculation
!
SUBROUTINE input_parameter_dyn()
  !
#if defined(MPI)
  USE mpi, only : MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_COMPLEX, &
  &               MPI_CHARACTER, MPI_LOGICAL
#endif
  USE shiftk_vals, ONLY : nomega, outrestart, calctype, z, e_max, e_min, stdout, myrank
  !
  IMPLICIT NONE
  !
  INTEGER :: iomega
#if defined(MPI)
  INTEGER ierr
#endif
  COMPLEX(8) :: omegamax, omegamin
  NAMELIST /dyn/ omegamax, omegamin, nomega, calctype, outrestart
  !
  calctype = "normal"
  nomega = 10
  omegamin = CMPLX(e_min, 0.01d0*(e_max - e_min), KIND(0d0))
  omegamax = CMPLX(e_max, 0.01d0*(e_max - e_min), KIND(0d0))
  outrestart = .FALSE.
  !
  IF(myrank == 0) READ(*,dyn,err=100)
  !
#if defined(MPI)
  call MPI_BCAST(nomega,   1, MPI_INTEGER,        0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(omegamin, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(omegamax, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(calctype, 256, MPI_CHARACTER,    0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(outrestart, 1, MPI_LOGICAL,      0, MPI_COMM_WORLD, ierr)
#endif
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Parameter for Spectrum  ##########"
  WRITE(stdout,*)
  WRITE(stdout,*) "           Max. of Omega : ", omegamax
  WRITE(stdout,*) "           Min. of Omega : ", omegamin
  WRITE(stdout,*) "           Num. of Omega : ", nomega
  WRITE(stdout,'(a,a)') "         Calculation type : ", calctype
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
#if defined(MPI)
  CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  STOP
  !
END SUBROUTINE input_parameter_dyn
!
! Input Hamiltonian with the Matrix Market Format
!
SUBROUTINE input_hamiltonian()
  !
#if defined(MPI)
  USE mpi, only : MPI_COMM_WORLD
#endif
  USE shiftk_vals, ONLY : ndim, inham, lBiCG, almost0, nproc, stdout
  USE ham_vals, ONLY : ham, ham_indx, nham, ndiag
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ndim2, iham, ham_indx0(2)
#if defined(MPI)
  INTEGER ierr
#endif
  REAL(8) :: ham_r, ham_i
  COMPLEX(8) :: ham0
  CHARACTER(100) :: ctmp
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Hamiltonian  ##########"
  WRITE(stdout,*)
  !
  IF(nproc /= 1) THEN
     WRITE(*,*) "ERROR ! MPI is not available in this mode."
#if defined(MPI)
     CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
     STOP
  END IF
  !
  OPEN(fi, file = TRIM(inham))
  READ(fi, *) ctmp
  READ(fi,*) ndim, ndim2, nham
  WRITE(stdout,*) "          Dim. of Hamiltonian : ", ndim, ndim2
  WRITE(stdout,*) "  Num. of Non-Zero Components : ", nham
  !
  IF(ndim2 /= ndim) THEN
     WRITE(stdout,*) "ERROR ! Hamiltonian is not square."
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
     WRITE(stdout,*) "  BiCG mathod is used."
     lBiCG = .TRUE.
  ELSE
     WRITE(stdout,*) "  COCG mathod is used."
     lBiCG = .FALSE.
  END IF
  !
END SUBROUTINE input_hamiltonian
!
! Input Right Hand Side Vector
!
SUBROUTINE input_rhs_vector()
  !
#if defined(MPI)
  USE mpi, only : MPI_COMM_WORLD
#endif
  USE shiftk_vals, ONLY : ndim, rhs, invec, e_min, e_max, nproc, stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ndim2, idim
#if defined(MPI)
  INTEGER ierr
#endif
  REAL(8) :: rhs_r, rhs_i
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Right Hand Side Vector ##########"
  WRITE(stdout,*)
  !
  IF(nproc /= 1) THEN
     WRITE(*,*) "ERROR ! MPI is not available in this mode."
#if defined(MPI)
     CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
     STOP
  END IF
  !
  OPEN(fi, file = TRIM(invec))
  READ(fi,*) ndim2
  WRITE(stdout,*) "  Dim. of RHS vector : ", ndim2
  !
  IF(ndim2 /= ndim) THEN
     WRITE(stdout,*) "ERROR ! Dimension is Incorrect."
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
  e_min = 0d0
  e_max = 1d0
  !
END SUBROUTINE input_rhs_vector
!
! Input restart variables from file
!
SUBROUTINE input_restart_parameter()
  !
#if defined(MPI)
  USE mpi, only : MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_COMPLEX
#endif
  USE shiftk_vals, ONLY : iter_old, alpha, beta, z_seed, nl, r_l_save, myrank, stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, iter, il
#if defined(MPI)
  INTEGER ierr
#endif
  REAL(8) :: z_seed_r, z_seed_i, alpha_r, alpha_i, beta_r, beta_i
  REAL(8) :: r_l_save_r, r_l_save_i
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Restart Parameter  ##########"
  WRITE(stdout,*)
  !
  IF(myrank == 0) THEN
     !
     OPEN(fi, file = 'output/TriDiagComp.dat')
     !
     READ(fi,*) iter_old
     WRITE(stdout,*) "  Num. of Iteration (Previous Run) : ", iter_old
     !
     ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl, iter_old))
     !
     READ(fi,*) z_seed_r, z_seed_i
     WRITE(stdout,*) "  Previous Omega_Seed : ", z_seed_r, z_seed_i
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
  END IF ! (myrank == 0)
  !
#if defined(MPI)
  call MPI_BCAST(iter_old, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  IF(myrank /= 0) THEN
     ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl, iter_old))
  END IF
  call MPI_BCAST(z_seed,   1,             MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(alpha,    iter_old,      MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(beta,     iter_old,      MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(r_l_save, iter_old * nl, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
#endif
  !
END SUBROUTINE input_restart_parameter
!
! Input restart variables from file
!
SUBROUTINE input_restart_vector()
  !
#if defined(MPI)
  USE mpi, only : MPI_COMM_WORLD
#endif
  USE shiftk_vals, ONLY : ndim, v2, v12, v4, v14, lBiCG, myrank, stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: fi = 10, ndim2, idim
  REAL(8) :: v2_r, v2_i, v12_r, v12_i
  CHARACTER(256) :: fname, cmyrank
#if defined(MPI)
  INTEGER ierr
#endif
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Input Restart Vector  ##########"
  WRITE(stdout,*)
  !
  WRITE(cmyrank,*) myrank
  WRITE(fname,'(a,a)') 'output/ResVec.dat', TRIM(ADJUSTL(cmyrank))
  OPEN(fi, file = TRIM(fname))
  !
  READ(fi,*) ndim2
  WRITE(stdout,*) "  Dim. of Residual vector : ", ndim2
  !
  IF(ndim2 /= ndim) THEN
     WRITE(stdout,*) "ERROR ! Dimension is Incorrect."
#if defined(MPI)
     CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
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
  USE shiftk_vals, ONLY : iter_old, alpha, beta, z_seed, r_l_save, nl, myrank, stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, iter, il
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Output Restart Parameter  ##########"
  WRITE(stdout,*)
  !
  IF(myrank == 0) THEN
     !
     OPEN(fo, file = 'output/TriDiagComp.dat')
     !
     WRITE(fo,*) iter_old
     WRITE(stdout,*) "  Num. of Iteration (Current Run) : ", iter_old
     !
     WRITE(stdout,'(a,2e13.5)') "   Current Omega_Seed : ", z_seed
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
  END IF !(myrank == 0)
  !
END SUBROUTINE output_restart_parameter
!
! Input restart variables from file
!
SUBROUTINE output_restart_vector()
  !
  USE shiftk_vals, ONLY : ndim, v2, v12, v4, v14, lBiCG, myrank, stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: fo = 20, idim
  CHARACTER(256) :: fname, cmyrank
  !
  WRITE(stdout,*)
  WRITE(stdout,*) "##########  Output Restart Vector  ##########"
  WRITE(stdout,*)
  !
  WRITE(cmyrank,*) myrank
  WRITE(fname,'(a,a)') 'output/ResVec.dat', TRIM(ADJUSTL(cmyrank))
  OPEN(fo, file = TRIM(fname))
  !
  WRITE(fo,*) ndim
  WRITE(stdout,*) "  Dim. of Residual vector : ", ndim
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
  USE shiftk_vals, ONLY : x_l, z, nomega, myrank
  !
  IMPLICIT NONE
  !
  INTEGER :: iz, fo = 20
  !
  IF(myrank == 0) THEN
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
  END IF ! (myrank == 0)
  !
END SUBROUTINE output_result
!
! Output Resulting Spectrum, True residual
!
SUBROUTINE output_result_debug()
  !
  USE shiftk_vals, ONLY : v2, ndim, x_l, rhs, z, nomega, myrank, stdout
  USE ham_prod_mod, ONLY : ham_prod
  USE lobpcg_mod, ONLY : zabsmax, zdotcMPI
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
  WRITE(stdout,*)
  WRITE(stdout,*) "#####  Check Results  #####"
  WRITE(stdout,*)
  !
  WRITE(stdout,*) "  Residual Vector"
  !
  DO iz = 1, nomega
     !
     CALL ham_prod(x_l(1:ndim,iz), v2)
     v2(1:ndim) = z(iz) * x_l(1:ndim,iz) - v2(1:ndim) - rhs(1:ndim)
     !
     !write(*,form) v2(1:ndim)
     write(*,'(a,i5,a,2e13.5,a,e13.5)') "DEBUG (", iz, "), omega = ", z(iz), &
     &            ", Res. = ", zabsmax(v2(1:ndim), ndim)
     !
  END DO
  !
  IF (myrank == 0) THEN
     !
     OPEN(fo, file = "output/dynamicalG.dat")
     !
     DO iz = 1, nomega
        !
        Gii = zdotcMPI(ndim, rhs,x_l(1:ndim,iz))
        write(fo,'(4e13.5)') DBLE(z(iz)), AIMAG(z(iz)), DBLE(Gii), AIMAG(Gii)
        !
     END DO
     !
     CLOSE(fo)
     !
  END IF ! (myrank == 0)
  !
END SUBROUTINE output_result_debug
!
END MODULE shiftk_io
