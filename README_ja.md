<a name= "japanese">
日本語 / [English](README.md#english)

<img src="doc/figs/komega.png" width="300">

# What ? 

Shifted-Krylov部分空間法に基づくソルバーライブラリと,
それを用いてHamiltonianと励起状態ベクトルから動的Green関数を計算するミニアプリである.

# Download

 * 最新版のパッケージはこちらから.
 
   https://github.com/issp-center-dev/Komega/releases/download/v1.0.0/komega-1.0.0.tar.gz

 * リリースノート
 
   https://github.com/issp-center-dev/Komega/releases/tag/v1.0.0
 * 過去のリリース
 
   https://github.com/issp-center-dev/Komega/releases

# Prerequisite

 * fortranコンパイラ
 * BLASライブラリ
 * LAPACKライブラリ(ミニアプリのみ使用)
 * MPIライブラリ(Optional)
 
# Documents

 * Manual for the Library
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/library/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/ja/_build/latex/komega.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/library/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/en/_build/latex/komega.pdf))
 * Manual for the sample program
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/software/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/ja/_build/latex/shiftk.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/software/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/en/_build/latex/shiftk.pdf))

# Directory Tree

 * `app/`: ミニアプリ関連のディレクトリ
   * `src/`: ミニアプリソースコードのディレクトリ
   * `sample/`: ミニアプリサンプル用ディレクトリ
     * `Shiftk.nb`: テスト用Mathematicaノートブック(開発者向け)
     * `denovo/`: ハミルトニアンや右辺ベクトルの入力ファイルを用いない場合の例
     * `from_file/`: ハミルトニアンや右辺ベクトルの入力ファイルを用いる場合の例
 * `doc/`: ドキュメント用ディレクトリ
   * `index.html` : ドキュメントのトップページ
   * `library/`: ライブラリのドキュメントのディレクトリ
   * `software/`: ミニアプリのドキュメントのディレクトリ
 * `configure`: ビルド環境指定ファイル
 * `src/`: ライブラリ本体のディレクトリ
 * `test/`: ライブラリのテスト用ディレクトリ

# Install

もっともシンプルには次のとおりである.

 * コンソールで `$ ./configure --prefix=**PREFIX**; make; make install` とタイプする.
 * `**PREFIX**`で指定したディレクトリに以下のものが作られる.
   * `**PREFIX**/lib/`内: 共有および静的ライブラリ
   * `**PREFIX**/include/komega.h`: C/C++用ヘッダーファイル
   * `**PREFIX**/bin/Shiftk.out`: ミニアプリ

詳しくはマニュアルを参照のこと.

# ミニアプリのテスト

 * `app/sample/denovo/`もしくは`app/sample/from_file/`ディレクトリに移動する.
 * `**PREFIX**/bin/ShiftK.out namelist.def`とやる.
 * `dynamicalG.dat`などが作られると成功.
 * `namelist.def`の書式はマニュアル参照

# ライブラリの使用方法

## プログラム内での各ルーチンの呼び出し方

### fortran/C/C++の場合

マニュアル参照

## ライブラリのリンク方法

``` bash
$ gfortran myprog.f90 -L **PREFIX**/src/shared -lshiftk -lblas -I **PREFIX**/src/shared
$ gcc myprog.c -L **PREFIX**/src/shared -lshiftk -lblas -I **PREFIX**/src/shared
```
など.

動的リンクを行ったファイルを実行するときには,
環境変数`LD_LIBRARY_PATH`に`**PREFIX**/lib`ディレクトリを追加しておく必要がある.

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
     
