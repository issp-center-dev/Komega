<a name= "japanese"> 日本語 /  [English](#english)

# What ? 

Shifted-Krylov部分空間法に基づくソルバーライブラリと,
それを用いてHamiltonianと励起状態ベクトルから動的Green関数を計算するミニアプリである.

# Prerequisite

 * fortran コンパイラ
 * BLASライブラリ
 
# Docments

 * Manual for the Library
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/library/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/ja/_build/latex/komega.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/library/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/en/_build/latex/komega.pdf))
 * Manual for the sample program
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/software/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/ja/_build/latex/komega.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/software/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/en/_build/latex/komega.pdf))

# Files in this package

 * `app/`: ミニアプリ関連のディレクトリ
   * `src/`: ミニアプリソースコードのディレクトリ
     * `dyn_mod.f90`: 動的Green関数計算用のサブルーチン群
     * `ham_prod.f90`: Hamiltonian-vector積のサブルーチン群
     * `lobpcg_mod.f90` : LOBPCG法のルーチン
     * `makefile`: ミニアプリのビルド用Makeファイル
     * `shiftk.f90`: メインプログラム
     * `shiftk_io.f90`: 入出力関連のサブルーチン群
     * `shiftk_vals.f90`: ミニアプリ内部共通変数モジュール
   * `sample/`: ミニアプリサンプル用ディレクトリ
     * `Shiftk.nb`: テスト用Mathematicaノートブック(開発者向け)
     * `denovo/`: ハミルトニアンや右辺ベクトルの入力ファイルを用いない場合の例
       * `namelist.def`: ミニアプリテスト用入力パラメーターファイル
     * `from_file/`: ハミルトニアンや右辺ベクトルの入力ファイルを用いる場合の例
       * `namelist.def`: ミニアプリテスト用入力パラメーターファイル
       * `zvo_Excited.dat`: テスト用励起ベクトルファイル(入力)
       * `zvo_Ham.dat`: テスト用Hamiltonianファイル(入力)
 * `doc/`: ドキュメント用ディレクトリ
   * `ShiftKSoft.pdf` : ミニアプリのマニュアル
   * `komega.pdf` : ライブラリのマニュアル
   * `library`: ライブラリのドキュメントのディレクトリ
     * `komega.tex`: ライブラリのマニュアルのソース
   * `software`: ミニアプリのドキュメントのディレクトリ
     * `KrylovSoft_ver.0.1.eps`: フロー図
     * `ShiftKSoft.tex`: ミニアプリのマニュアルのソース
 * `make.sys`: ビルド環境指定ファイル
 * `makefile`: Makeファイル
 * `src/`: ライブラリのソースコードのディレクトリ
   * `makefile`: ライブラリのビルド用Makeファイル
   * `komega_bicg.f90`: Shifted BiCG法ライブラリ用サブルーチン群
   * `komega_cg_c.f90`: Shifted CG法(複素Hamiltonian)ライブラリ用サブルーチン群
   * `komega_cg_r.f90`: Shifted CG法(実Hamiltonian)ライブラリ用サブルーチン群
   * `komega_cocg.f90`: Shifted COCG法ライブラリ用サブルーチン群
   * `komega_math.f90`: BLASインターフェイスモジュール
   * `komega_vals.f90`: ライブラリ内部共通変数モジュール
 * `test/`: ライブラリのテスト用ディレクトリ
   * `krylov.in`: テスト用入力パラメーターファイル
   * `make_ham.f90`: 擬似Hamiltonianを乱数で生成するサブルーチン群
   * `makefile`: テスト用Makeファイル
   * `solve_cc.f90`: 複素Hamiltonian-複素振動数(BiCG)テスト
   * `solve_cr.f90`: 複素Hamiltonian-実振動数(CG-C)テスト
   * `solve_rc.f90`: 実Hamiltonian-複素振動数(COCG)テスト
   * `solve_rr.f90`: 実Hamiltonian-実振動数(CG-R)テスト

# Build

 * 必要に応じて`make.sys`を編集する.
   * `F90`:fortranコンパイルコマンド
   * `FFLAGS`: オプション.BLASのリンク(`-lblas`や`-mkl`)など.
   * `MPIF90`:MPI版fortranコンパイルコマンド. MPI版が必要ない場合には指定しなくて良い.
 * `$ make`

以下のものが作られる.

 * `src/libkomega.a` : シリアル版静的ライブラリ
 * `src/shared/libkomega.so` : シリアル版動的ライブラリ
 * `src/mpi/libpkomega.a` : MPI版静的ライブラリ(Optional)
 * `src/shared_mpi/libpkomega.so` : MPI版動的ライブラリ(Optional)
 * `app/src/Shiftk.out` : ミニアプリ
 * `app/src/mpi/Shiftk.out` : MPI版ミニアプリ (Oprional)

# ミニアプリのテスト

 * `app/sample/denovo/`もしくは`app/sample/from_file/`ディレクトリに移動する.
 * `$ ../../src/ShiftK.out namelist.def`とやる.
 * `dynamicalG.dat`などが作られると成功.
 * `namelist.def`の書式はマニュアル参照

# ライブラリの使用方法

## プログラム内での各ルーチンの呼び出し方

### fortran/C/C++の場合

マニュアル参照(作成中)

## ライブラリのリンク方法

### fortranの場合

静的リンク
```
$ ifort myprog.f90 -L パス/src -lshiftk -I パス/src
```
動的リンク
```
$ ifort myprog.f90 -L パス/src/shared -lshiftk -I パス/src/shared
```
動的リンクを行ったファイルを実行するときには,
環境変数`LD_LIBRARY_PATH`に`src/shared`ディレクトリを追加しておく必要がある.

### C/C++の場合

# ライブラリのテスト(Optional)

 * `test/`ディレクトリに移動する.
 * `$ ./solve_cc.x < krylov.in`とやる.
 * Residual vectorが十分小さくなっていれば成功
 * `solve_rc.x`なども同様
 * `krylov.in`(名称は自由)のパラメーターは次の通り
   * `ndim`: 擬似ハミルトニアンの次元
   * `nl`: 射影のテスト用. 解ベクトルの先頭から`nl(<=ndim)`番目までを計算する.
   * `nz`: 振動数の点数
   * `itermax`: 最大反復回数
   * `threshold`: 収束判定のthreshold
   * `rnd_seed`: 擬似Hamiltonian生成のための乱数の種
   * このnamelistの後ろに振動数を1行に1個ずつ書く.
     
----
   
[日本語]( #japanese ) / <a name= "english"> English 
  
# What ? 

This package provides the solver library based on Shifted-Krylov subspace method and the software to calculate dynamical Green function by inputting the Hamiltonian and the vector of the excited state.


# Prerequisite

 * fortran compiler
 * BLAS library  

# Docments

 * Manual for the Library
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/library/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/ja/_build/latex/komega.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/library/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/en/_build/latex/komega.pdf))
 * Manual for the sample program
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/software/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/ja/_build/latex/komega.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/software/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/en/_build/latex/komega.pdf))

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
  

