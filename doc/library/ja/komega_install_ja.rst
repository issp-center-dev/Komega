インストール方法
================

大まかな手順
------------

最もシンプルには次のとおりである.

.. code-block:: bash

   $ ./configure --prefix=**PREFIX**

これにより, ビルドに必要なコンパイラやライブラリ等の環境のチェックが行われ,
Makefile等が作成される.
``**PREFIX**`` はライブラリ等がインストールされるディレクトリを示す.
なにも指定しないと ``/use/local/`` が設定され, 後述の ``make install`` で
``/usr/local/lib`` 内にライブラリが置かれる (したがって, 管理者権限がない場合には ``**PREFIX**`` を
別の場所に指定しなければならない).
``configure`` にはこの他にも様々なオプションがあり,必要に応じて用途や環境に合わせてそれらを使用する.
詳しくは :ref:`configoption` を参照.

``configure`` の実行が正常に行われ, ``Makefile`` が生成された後は

.. code-block:: bash

   $ make

とタイプしてライブラリ等のビルドを行う.これが成功したのちに

.. code-block:: bash

   $ make install

とすると, ライブラリが ``**PREFIX**/lib`` に, ミニアプリが ``**PREFIX**/bin`` に置かれる.
``make install`` をしなくても, ビルドをしたディレクトリ内にあるライブラリやミニアプリを使うことは可能であるが,
使い勝手がやや異なる.

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

``--disable-zdot``

   デフォルト: ``--enable-zdot``.
   MacOSXの標準のBLAS等では, ZDOTCおよびZDOTU関数が正常に動作しないため,
   このオプションでこれらの関数を使わないようにする.

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
