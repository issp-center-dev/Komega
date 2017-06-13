<a name= "english">
[日本語]( README_ja.md#japanese ) / English 

<img src="doc/figs/komega.png" width="300">

# What ? 

This package provides the solver library based on Shifted-Krylov subspace method and the software to calculate dynamical Green function by inputting the Hamiltonian and the vector of the excited state.

# Download

 * Latest release

   https://github.com/issp-center-dev/Komega/releases/download/v2.0.0/komega-2.0.0.tar.gz
 * Release note

   https://github.com/issp-center-dev/Komega/releases/tag/v2.0.0
 * Other version
 
   https://github.com/issp-center-dev/Komega/releases
   
# Prerequisite

 * fortran compiler
 * BLAS library  
 * LAPACK library (Used in the sample program)
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
 * `configure`: Configuration script to build
 * `src/`: The source directory for the libraries
 * `test/`: The test directory for the libraries

# Install

The simplest procedure is as follows:

 * Type `$ ./configure --prefix=install_dir; make; make install`
   where `install_dir` should be replaced with the full path of the directory where
   the library will be stored.
 * The following objects are generated in the directory specified by `install_dir`.
   * In `install_dir/lib/`: Static and shared libraries.
   * In `install_dir/include/`: Header file for C/C++.
   * `install_dir/bin/Shiftk.out`: Sample program

For more details, please see the manual.

# Test of the software

 * Change the directory to `app/sample/denovo/` or `app/sample/from_file/`.
 * Type the command `$ install_dir/bin/ShiftK.out namelist.def`.
 * When the software works well, the files such as `dynamicalG.dat` will be generated.
 * The details of the file format of `namelist.def` is written in the manual.

# Usage of libraries

## How to call each routines in the program

### For fortran/C/C++

See the manual.

## How to link the libraries

```
$ gfortran myprog.f90 -L install_dir/lib -lshiftk -lblas -I install_dir/include
$ gcc myprog.c -L install_dir/lib -lshiftk -lblas -I install_dir/include
```
etc.

Add the `install_dir/lib` directry to the environment
variable `LD_LIBRARY_PATH` to execute the file with dynamic link.

# Test for libraries(Optional)

 * Move to the `test/` directory.
 * Type the command `$ ./solve_cc.x < complex_freq.in`.
 * The program is normally finished when Residual vector becomes sufficiently small. 
 * You can check other programs (`solve_rc.x`, `solve_cr.x`, `solve_rr.x`) in a similar way.
   For `solve_cr.x` and `solve_rr.x, please use `real_freq.in` as an input file.
 * The parameters in `krylov.in`(you can modify the file name freely) are as follows
   * `ndim`: The dimeinsion of the psuedo Hamiltonian
   * `nl`: This parameter is used to test the projection.
     The vectors are calculated from the target vector up to `nl(<=ndim)`th vector.
   * `nz`: The number of the frequences to calculate.
   * `itermax`: The maximum number of iterations.
   * `threshold`: The threshold to judge the convergence.
   * `rnd_seed`: The seed of random number to generate pseudo Hamiltonian.
   * Write the value of each frequencies line by line after this namelist.
