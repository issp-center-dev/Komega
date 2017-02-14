<a name= "english">
[日本語]( README_ja.md#japanese ) / English 

<img src="doc/figs/komega.png" width="300">

# What ? 

This package provides the solver library based on Shifted-Krylov subspace method and the software to calculate dynamical Green function by inputting the Hamiltonian and the vector of the excited state.

# Download

 * Latest release

   https://github.com/issp-center-dev/Komega/releases/download/v1.0.0/komega-1.0.0.tar.gz
 * Release note

   https://github.com/issp-center-dev/Komega/releases/tag/v1.0.0
 * Other version
 
   https://github.com/issp-center-dev/Komega/releases
   
# Prerequisite

 * fortran compiler
 * BLAS library  
 * LAPCK library  
 * MPI library (Optional)

# Docments

 * Manual for the Library
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/library/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/ja/_build/latex/komega.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/library/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/en/_build/latex/komega.pdf))
 * Manual for the sample program
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/software/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/ja/_build/latex/shiftk.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/software/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/en/_build/latex/shiftk.pdf))

# Directory Tree

 * `app/`: The directory for the software
   * `src/`: The source directory for the software
   * `sample/`: The sample directory
     * `Shiftk.nb`: The mathmatica note book to test the software (for developers)
     * `denovo/`: Sample files not to input either Hamiltonian and the right-hand side initial vector
     * `from_file/`: Sample files to input both Hamiltonian and the right-hand side initial vector
 * `doc/`: The documents directory (only japanese, english version will be provided from ver.1.0)
   * `index.html` : Top page for documents
   * `library/`: Directory for the document of the librariy
   * `software/`: Directory for the document of the sample program
 * `make.sys`: Configuration file to build
 * `makefile`: Makefile
 * `src/`: The source directory for the libraries
   * `mpi/`: Directory for the library with MPI
   * `shared/`: Directory for the dynamic library
   * `shared_mpi/`: Directory for the dynamic library with MPI
 * `test/`: The test directory for the libraries

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

- Static link
  ```
  $ gfortran myprog.f90 -L Path/src -lshiftk -I Path/src
  $ gcc myprog.c -L Path/src -lshiftk -I Path/src
  ```
  etc.

- Dynamical link
  ```
  $ gfortran myprog.f90 -L Path/src/shared -lshiftk -I Path/src/shared
  $ gcc myprog.c -L Path/src/shared -lshiftk -I Path/src/shared
  ```
  etc.

  Add the `src/shared` directry to the environment
  variable `LD_LIBRARY_PATH` to execute the file with dyncamic link.

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
  

