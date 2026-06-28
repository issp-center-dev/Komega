インストール方法
================

大まかな手順
------------

最もシンプルには次のとおりである.

.. code-block:: bash

   $ ./configure --prefix=install_dir

これにより, ビルドに必要なコンパイラやライブラリ等の環境のチェックが行われ,
Makefile等が作成される.
ただし ``install_dir`` はインストール先のディレクトリの絶対パスとする (以後各自のディレクトリ名で読み替えること).
なにも指定しないと ``/use/local/`` が設定され, 後述の ``make install`` で
``/usr/local/lib`` 内にライブラリが置かれる (したがって, 管理者権限がない場合には ``install_dir`` を
別の場所に指定しなければならない).
``configure`` にはこの他にも様々なオプションがあり,必要に応じて用途や環境に合わせてそれらを使用する.
詳しくは :ref:`configoption` を参照.

``configure`` の実行が正常に行われ, ``Makefile`` が生成された後は

.. code-block:: bash

   $ make

とタイプしてライブラリ等のビルドを行う.

必要に応じて, 次のコマンドで回帰テストを実行できる.

.. code-block:: bash

   $ make check

これはサンプルソルバーを実行し, 許容誤差内に収束することを確認する
(失敗した場合は非ゼロのステータスを返す).

ビルドが成功したのちに

.. code-block:: bash

   $ make install

とすると, ライブラリが ``install_dir/lib`` に, ミニアプリが ``install_dir/bin`` に置かれる.
``make install`` をしなくても, ビルドをしたディレクトリ内にあるライブラリやミニアプリを使うことは可能であるが,
使い勝手がやや異なる.

共有リンクを行ったプログラムの実行時にライブラリを探しにいけるよう,
環境変数 ``LD_LIBRARY_PATH`` に :math:`K\omega` をインストールしたディレクトリを追加する.

.. code-block:: bash

   $ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:install_dir/lib

.. _configoption:

configureのオプション
---------------------

configureには多数のオプションと変数があり, それらを組み合わせて指定する.
指定しない場合にはデフォルト値が使われる.

.. code-block:: bash

  $ ./configure --prefix=/home/komega/ --with-mpi=yes FC=mpif90

おもなものを次に挙げる.

``---prefix``

   デフォルト: ``---prefix=/usr/local/``.
   ライブラリ等のインストールを行うディレクトリツリーを指定する.

``--with-mpi``

   デフォルト: ``--with-mpi=no`` (MPIを用いない).
   MPIを用いるか (``--with-mpi=yes``), 否かを指定する.

``--with-openmp``

   デフォルト: ``--with-openmp=yes`` (OpenMPを用いる).
   OpenMPを用いるか否か (``--with-openmp=no``) を指定する.

``--enable-shared``

   デフォルト: ``--enable-shared``.
   共有ライブラリを作成するか否か

``--enable-static``

   デフォルト: ``--enable-static``.
   静的ライブラリを作成するか否か.

``--enable-zdot`` / ``--disable-zdot``

   デフォルト: 自動判定.
   BLASの ``ZDOTC`` / ``ZDOTU`` 関数を使用するかどうかを制御する.
   一部のBLAS (特にmacOSの標準のBLASやAccelerateフレームワーク) では,
   これらの関数が複素数の結果を :math:`K\omega` が期待するレジスタ渡しではなく,
   隠れた第一引数を介して返す (``f2c``/``g77`` の呼出規約) ため,
   プログラムがクラッシュする.
   デフォルトでは ``configure`` が小さなテストプログラムをコンパイル・実行して
   これを自動判定する. BLASの ``ZDOTC`` が正しく動作する場合はそれを使用し,
   そうでない場合はFortran組込み関数 ``DOT_PRODUCT`` / ``SUM`` に切り替える
   (マクロ ``-D__NO_ZDOT`` を定義するのと等価).
   ``--enable-zdot`` を指定するとBLASの関数の使用を強制し,
   ``--disable-zdot`` を指定すると組込み関数へのフォールバックを強制する
   (いずれも自動判定を行わない).
   クロスコンパイル時 (テストプログラムを実行できない場合) には,
   安全側のフォールバック (``--disable-zdot`` 相当) が自動的に選択される.

``--enable-threadsafe``

   デフォルト: ``--disable-threadsafe``.
   もしもOpenMPのパラレルリージョンの内側で :math:`K\omega` を呼び出したい
   (それぞれのスレッドで異なる問題を解きたい)場合には,
   このオプションを用いる (**試験的**).

``FC``

   デフォルト: システムにインストールされているfortranコンパイラをスキャンして,
   自動的に設定する. ``--with-mpiP`` を指定した時にはそれに応じたコマンド
   (``mpif90`` 等)を自動で探し出して設定する. 
   ``configure`` の最後に出力される ``FC`` が望んだものでは無かった場合には
   ``./configure FC=gfortran`` のように手で指定する.

``--help``

   このオプションを指定した時には, ビルドの環境設定は行われず,
   上記を含めたすべてのオプションを表示する.
