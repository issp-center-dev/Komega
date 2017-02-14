<a name= "english">
[日本語]( README_ja.md#japanese ) / English 
  
# What ? 

This package provides the solver library based on Shifted-Krylov subspace method and the software to calculate dynamical Green function by inputting the Hamiltonian and the vector of the excited state.

# [Download](https://github.com/issp-center-dev/Komega/releases)

# Prerequisite

 * fortran compiler
 * BLAS library  

# Docments

 * Manual for the Library
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/library/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/ja/_build/latex/komega.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/library/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/en/_build/latex/komega.pdf))
 * Manual for the sample program
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/software/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/ja/_build/latex/shiftk.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/software/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/en/_build/latex/shiftk.pdf))

# Files in this package

 * `app/`: The directory for the software
   * `src/`: The source directory for the software
     * `dyn_mod.f90`: The subroutines to calculate the dynamical Green's functions
     * `ham_prod.f90`: Subroutines to mltply the Hamiltonian by the vector
     * `lobpcg_mod.f90` : Routines for the LOBPCG method
     * `makefile`: Makefile to bulid the software
     * `shiftk.f90`: Main program
     * `shiftk_io.f90`: Subroutines for the input/output data
     * `shiftk_vals.f90`: Modules for the internal common variables
   * `sample/`: The sample directory
     * `Shiftk.nb`: The mathmatica note book to test the software (for developers)
     * `denovo/`: Sample files not to input either Hamiltonian and the right-hand side initial vector
       * `namelist.def`: The input parameter file to test the software
     * `from_file/`: Sample files to input both Hamiltonian and the right-hand side initial vector
       * `namelist.def`: The input parameter file to test the software
       * `zvo_Excited.dat`: The input file of the initial excited vector to test the software
       * `zvo_Ham.dat`: The input file of the Hamiltonian to test the software
 * `doc/`: The documents directory (only japanese, english version will be provided from ver.1.0)
   * `ShiftKSoft.pdf` : The manual for the software 
   * `komega.pdf` : The manual for the libraries
   * `library`: The documents directory for the libraries
     * `komega.tex`: The source file of the manual for the libraries 
   * `software`: the documents directory for the software
     * `KrylovSoft_ver.0.1.eps`: Flow diagram
     * `ShiftKSoft.tex`: The source file of the manual for the software
 * `make.sys`: Configuration file to build
 * `makefile`: Makefile
 * `src/`: The source directory for the libraries
   * `makefile`: Makefile to build the libraries
   * `komega_bicg.f90`: Subroutines of the libraries for Shifted BiCG method
   * `komega_cg_c.f90`: Subroutines of the libraries for Shifted CG method (complex Hamiltonian)
   * `komega_cg_r.f90`: Subroutines of the libraries for Shifted CG method (real Hamiltonian)
   * `komega_cocg.f90`: Subroutines of the libraries for Shifted COCG method
   * `komega_math.f90`: Interface modules of BLAS
   * `komega_vals.f90`: Modules for the internal common variables
 * `test/`: The test directory for the libraries
   * `krylov.in`: The initial parameter files to test the libraries
   * `make_ham.f90`: Subroutines to generate the sample Hamiltonian by using the random number
   * `makefile`: Makefile to test the libraries
   * `solve_cc.f90`: The test file in the case of the complex Hamiltonian and the complex frequencies (BiCG method)
   * `solve_cr.f90`: The test file in the case of the complex Hamiltonian and the real frequencies (CG-C method)
   * `solve_rc.f90`: The test file in the case of the real Hamiltonian and the complex frequencies (COCG method)
   * `solve_rr.f90`: The test file in the case of the real Hamiltonian and the real frequencies (CG-R method)

# Build

 * Edit `make.sys` for your circumstances.
   * `F90`:The fortran compile command
   * `FFLAGS`:Options such as the BLAS link; `-lblas`, `-mkl`.
   * `MPIF90`:The MPI fortran compile command. If you do not use the MPI, you can ignore this command.
 * `$ make`

The following objects are generated.

 * `src/libkomega.a` : Static library in serial version
 * `src/shared/libkomega.so` : Dynamical library in serial version
 * `src/mpi/libpkomega.a` : Static library (MPI) in serial version (Optional)
 * `src/shared_mpi/libpkomega.so` : Dynamical library (MPI) in serial version (Optional)
 * `app/src/Shiftk.out` : Software
 * `app/src/mpi/Shiftk.out` : Software in MPI version(Optional)

# Test of the software

 * Change the directory to `app/sample/denovo/` or `app/sample/from_file/`.
 * Type the command `$ ../../src/ShiftK.out namelist.def`.
 * When the software works well, the files such as `dynamicalG.dat` will be generated.
 * The details of the file format of `namelist.def` is written in the manual.

# Usage of libraries

## How to call each rouchines in the program

### For fortran/C/C++

See the manual (the explanation will be added in ver.0.1)

## How to link the libraries

### For fortran

- Static link
```
$ ifort myprog.f90 -L Path/src -lshiftk -I Path/src
```

- Dynamical link
```
$ ifort myprog.f90 -L Path/src/shared -lshiftk -I Path/src/shared
```  
Add the `src/shared` directry to the environment variable `LD_LIBRARY_PATH` to execute the file with dyncamic link.


### For C/C++

# Test for libraries(Optional)

 * Move to the `test/` directory.
 * Type the command `$ ./solve_cc.x < krylov.in`.
 * The program is normally finished when Residual vector becomes sufficiently small. 
 * You can check other programs such as `solve_rc.x` in a similar way.
 * The parameters in `krylov.in`(you can modify the file name freely) are as follows
   * `ndim`: The dimeinsion of the psuedo Hamiltonian
   * `nl`: This parameter is used to test the projection. The vectors are calculated from the target vector up to `nl(<=ndim)`th vector.
   * `nz`: The number of the frequences to calculate.
   * `itermax`: The maximum number of iterations.
   * `threshold`: The threshold to judge the convergence.
   * `rnd_seed`: The seed of random number to generate pseudo Hamiltonian.
   * Write the value of each frequencies line by line after this namelist.
  

