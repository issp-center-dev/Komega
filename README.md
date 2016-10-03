# What ?

Shifted-Krylov部分空間法に基づくソルバーライブラリと,
それを用いてHamiltonianと励起状態ベクトルから動的Green関数を計算するミニアプリである.

# Prerequisite

 * fortran コンパイラ
 * BLASライブラリ

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
   * `library`: ライブラリのドキュメントのディレクトリ
     * `ShiftK.tex`: ライブラリの設計書兼マニュアル
   * `software`: ミニアプリのドキュメントのディレクトリ
     * `KrylovSoft_ver.0.1.eps`: フロー図
     * `ShiftKSoft.tex`: ミニアプリの設計書兼マニュアル
 * `make.sys`: ビルド環境指定ファイル
 * `makefile`: Makeファイル
 * `src/`: ライブラリのソースコードのディレクトリ
   * `makefile`: ライブラリのビルド用Makeファイル
   * `shifted_bicg.f90`: Shifted BiCG法ライブラリ用サブルーチン群
   * `shifted_cg_c.f90`: Shifted CG法(実Hamiltonian)ライブラリ用サブルーチン群
   * `shifted_cg_r.f90`: Shifted CG法(複素Hamiltonian)ライブラリ用サブルーチン群
   * `shifted_cocg.f90`: Shifted COCG法ライブラリ用サブルーチン群
   * `shifted_krylov_math.f90`: BLASインターフェイスモジュール
   * `shifted_krylov_vals.f90`: ライブラリ内部共通変数モジュール
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

 * `src/libshiftk.a` : シリアル版静的ライブラリ
 * `src/shared/libshiftk.so` : シリアル版動的ライブラリ
 * `src/mpi/libpshiftk.a` : MPI版静的ライブラリ(Optional)
 * `src/shared_mpi/libpshiftk.so` : MPI版動的ライブラリ(Optional)
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
