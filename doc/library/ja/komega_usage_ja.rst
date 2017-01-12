使用方法
========

各ライブラリともユーザーはライブラリ名および型を指定し,

-  初期設定 (:ref:`init`)

-  アップデート (:ref:`update`)

-  (オプション) 再計算用の情報を取り出す. (:ref:`getcoef`, :ref:`getvec`)

-  終了関数 (:ref:`finalize`)

の手順で関数を使用することで, 計算が実施される. なお,
リスタートを行う場合には

-  前回の計算で残した再計算用の情報を用いた初期設定(:ref:`restart`)

-  アップデート (:ref:`update`)

-  (オプション) 更なる再計算用の情報を取り出す. (:ref:`getcoef`, :ref:`getvec`)

-  終了関数 (:ref:`finalize`)

の手順で実行する. fortran から呼び出すときには

.. code-block:: fortran

   USE komega_cg_r ! 実ベクトルに対する共役勾配法
   USE komega_cg_c ! 複素ベクトルに対する共役勾配法
   USE komega_cocg ! 共線直交共役勾配法
   USE komega_bicg ! 双共役勾配法

のようにモジュールを呼び出す(すべてのモジュールを呼び出す必要はなく,
行う計算の種類に対応するものだけでよい). 
MPI/Hybrid並列版のルーチンを利用するときには,

.. code-block:: fortran

   USE pkomega_cg_r
   USE pkomega_cg_c
   USE pkomega_cocg
   USE pkomega_bicg

のようにする.

C/C++で書かれたプログラムから呼び出すときには、

.. code-block:: c

    #include komega_cg_r.h
    #include komega_cg_c.h
    #include komega_cocg.h
    #include komega_bicg.h

のようにヘッダーファイルを読み込む。
また、スカラー引数はすべてポインタとして渡す。
MPI/Hybrid並列版のルーチンを利用するときには,

.. code-block:: c

    #include pkomega_cg_r.h
    #include pkomega_cg_c.h
    #include pkomega_cocg.h
    #include pkomega_bicg.h

のようにする。
またライブラリに渡すコミュニケーター変数を、次のようにC/C++のものからfortranのものに変換する。

.. code-block:: c

      comm_f = MPI_Comm_c2f(comm_c);

各ルーチンの説明
----------------

.. _init:

\*_init
~~~~~~~

ライブラリ内部変数の割り付けおよび初期化を行う.
シフト線形方程式を解く前に, 一番初めに実行する.

構文

   Fortran シリアル/OpenMP版

   .. code-block:: fortran

      CALL komega_cg_r_init(ndim, nl, nz, x, z, itermax, threshold)
      CALL komega_cg_c_init(ndim, nl, nz, x, z, itermax, threshold)
      CALL komega_cocg_init(ndim, nl, nz, x, z, itermax, threshold)
      CALL komega_bicg_init(ndim, nl, nz, x, z, itermax, threshold)

   Fortran MPI/Hybrid並列版

   .. code-block:: fortran

      CALL pkomega_cg_r_init(ndim, nl, nz, x, z, itermax, threshold, comm)
      CALL pkomega_cg_c_init(ndim, nl, nz, x, z, itermax, threshold, comm)
      CALL pkomega_cocg_init(ndim, nl, nz, x, z, itermax, threshold, comm)
      CALL pkomega_bicg_init(ndim, nl, nz, x, z, itermax, threshold, comm)

   C/C++ シリアル/OpenMP版

   .. code-block:: c

      komega_cg_r_init(&ndim, &nl, &nz, x, z, &itermax, &threshold);
      komega_cg_c_init(&ndim, &nl, &nz, x, z, &itermax, &threshold);
      komega_cocg_init(&ndim, &nl, &nz, x, z, &itermax, &threshold);
      komega_bicg_init(&ndim, &nl, &nz, x, z, &itermax, &threshold);

   C/C++ MPI/Hybrid並列版

   .. code-block:: c

      pkomega_cg_r_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);
      pkomega_cg_c_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);
      pkomega_cocg_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);
      pkomega_bicg_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);

パラメーター

   .. code-block:: fortran

      INTEGER,INTENT(IN) :: ndim
   ..

      線形方程式の次元.
      以降のサブルーチンのパラメーターの次元で現れる ``ndim`` は
      これと同じものになる.

   .. code-block:: fortran

      INTEGER,INTENT(IN) :: nl
   ..

      射影された解ベクトルの次元.
      以降のサブルーチンのパラメーターの次元で現れる ``nl`` は
      これと同じものになる.

   .. code-block:: fortran
                
      INTEGER,INTENT(IN) :: nz
   ..

      シフト点の数.
      以降のサブルーチンのパラメーターの次元で現れる ``nz`` は
      これと同じものになる.

   .. code-block:: fortran

      REAL(8),INTENT(OUT) :: x(nl*nz) ! ("CG_R_init", "cg_c_init" の場合)
      COMPLEX(8),INTENT(OUT) :: x(nl*nz) ! (それ以外)
   ..

      解ベクトル. ``0`` ベクトルが返される.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: z(nz) ! ("CG_R_init", "cg_c_init" の場合)
      COMPLEX(8),INTENT(IN) :: z(nz) ! (それ以外)
   ..

      シフト点.
         
   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: itermax
   ..

      リスタート用配列の割り付けのための最大反復回数.
      これを ``0`` にした場合にはリスタート用配列を割りつけない
      (したがって後述のリスタート用変数の出力を行えない)

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: threshold
   ..

      収束判定用しきい値.
      シード方程式の残差ベクトルの2-ノルムがこの値を下回った時に収束したと判定する.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: comm
   ..

      MPI/Hybrid並列版のみ.
      MPIのコミニュケーター( ``MPI_COMM_WORLD`` など)を入れる.

.. _restart:
   
\*_restart
~~~~~~~~~~

リスタートを行う場合に :ref:`init` の代わりに用いる.
ライブラリ内部変数の割り付けおよび初期化を行う.
シフト線形方程式を解く前に, 一番初めに実行する.

構文

   Fortran (シリアル/OpenMP版)

   .. code-block:: fortran

      CALL komega_cg_r_restart(ndim, nl, nz, x, z, itermax, threshold, status, &
      &                 iter_old, v2, v12, alpha_save, beta_save, z_seed, r_l_save)
      CALL komega_cg_c_restart(ndim, nl, nz, x, z, itermax, threshold, status, &
      &                 iter_old, v2, v12, alpha_save, beta_save, z_seed, r_l_save)
      CALL komega_cocg_restart(ndim, nl, nz, x, z, itermax, threshold, status, &
      &                 iter_old, v2, v12, alpha_save, beta_save, z_seed, r_l_save)
      CALL komega_bicg_restart(ndim, nl, nz, x, z, itermax, threshold, status, &
      &                 iter_old, v2, v12, v4, v14, alpha_save, beta_save, &
      &                 z_seed, r_l_save)

   Fortran (MPI/ハイブリッド並列版)

   .. code-block:: fortran
                   
      CALL pkomega_cg_r_restart(ndim, nl, nz, x, z, itermax, threshold, comm, status, &
      &                 iter_old, v2, v12, alpha_save, beta_save, z_seed, r_l_save)
      CALL pkomega_cg_c_restart(ndim, nl, nz, x, z, itermax, threshold, comm, status, &
      &                 iter_old, v2, v12, alpha_save, beta_save, z_seed, r_l_save)
      CALL pkomega_cocg_restart(ndim, nl, nz, x, z, itermax, threshold, comm, status, &
      &                 iter_old, v2, v12, alpha_save, beta_save, z_seed, r_l_save)
      CALL pkomega_bicg_restart(ndim, nl, nz, x, z, itermax, threshold, comm, status, &
      &                 iter_old, v2, v12, v4, v14, alpha_save, beta_save, &
      &                 z_seed, r_l_save)

   C/C++ (シリアル/OpenMP版)

   .. code-block:: c

      komega_cg_r_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, status, &
      &                 &iter_old, v2, v12, alpha_save, beta_save, &z_seed, r_l_save);
      komega_cg_c_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, status, &
      &                 &iter_old, v2, v12, alpha_save, beta_save, &z_seed, r_l_save);
      komega_cocg_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, status, &
      &                 &iter_old, v2, v12, alpha_save, beta_save, &z_seed, r_l_save);
      komega_bicg_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, status, &
      &                 &iter_old, v2, v12, v4, v14, alpha_save, beta_save, &
      &                 &z_seed, r_l_save);

   C/C++ (MPI/ハイブリッド並列版)

   .. code-block:: c

      pkomega_cg_r_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm, status, &
      &                 &iter_old, v2, v12, alpha_save, beta_save, &z_seed, r_l_save);
      pkomega_cg_c_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm, status, &
      &                 &iter_old, v2, v12, alpha_save, beta_save, &z_seed, r_l_save);
      pkomega_cocg_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm, status, &
      &                 &iter_old, v2, v12, alpha_save, beta_save, &z_seed, r_l_save);
      pkomega_bicg_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm, status, &
      &                 &iter_old, v2, v12, v4, v14, alpha_save, beta_save, &
      &                 &z_seed, r_l_save);

パラメーター

   .. code-block:: fortran

      INTEGER,INTENT(IN) :: ndim
      INTEGER,INTENT(IN) :: nl
      INTEGER,INTENT(IN) :: nz
      REAL(8),INTENT(OUT) :: x(nl*nz)
      REAL(8),INTENT(IN) :: z(nz) ! ("CG_R_restart", "cg_c_restart" の場合)
      COMPLEX(8),INTENT(IN) :: z(nz) ! (それ以外)
      INTEGER,INTENT(IN) :: itermax
      REAL(8),INTENT(IN) :: threshold
      INTEGER,INTENT(IN) :: comm
   ..
   
      :ref:`init` と同様.

   .. code-block:: fortran

      INTEGER,INTENT(OUT) :: status(3)
   ..
   
      エラーコードを返す.

      第一成分( ``status(1)``)
      
         解が収束した場合,
         もしくは計算が破綻した場合には現在の総反復回数に
         マイナスが付いた値が返される.
         それ以外の場合には現在の総反復回数(マイナスが付かない)が返される.
         ``status(1)`` が正の値の時のみ反復を続行できる.
         それ以外の場合は反復を進めても有意な結果は得られない.

      第二成分( ``status(2)``)
   
         ``itermax`` を有限にして, かつ ``itermax`` 回の反復で
         収束に達しなかった場合には ``1`` が返される.
         :math:`\alpha` が発散した場合には ``2`` が返される.
         :math:`\pi_{\rm seed}` が0にになった場合には ``3`` が返される.
         ``COCG_restart`` もしくは ``BiCG_restart`` で,
         残差ベクトルと影の残差ベクトルが直交した場合には ``4`` が返される.
         それ以外の場合には ``0`` が返される.

      第三成分( ``status(3)``)
      
         シード点のindexが返される.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: iter_old
   ..
   
      先行する計算での反復回数.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: v2(ndim) ! ("CG_R_restart" の場合)
      COMPLEX(8),INTENT(IN) :: v2(ndim) ! (それ以外)
   ..
   
      先行する計算での最後の残差ベクトル.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: v12(ndim) ! ("CG_R_restart" の場合)
      COMPLEX(8),INTENT(IN) :: v12(ndim) ! (それ以外)
   ..

      先行する計算での最後から2番目の残差ベクトル.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: alpha_save(iter_old) ! ("CG_R_restart", "cg_c_restart"の場合)
      COMPLEX(8),INTENT(IN) :: alpha_save(iter_old) ! (それ以外)
   ..                   

      先行する計算での各反復での(Bi)CG法のパラメーター :math:`\alpha`.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: beta_save(iter_old) ! ("CG_R_restart", "cg_c_restart"の場合)
      COMPLEX(8),INTENT(IN) :: beta_save(iter_old) ! (それ以外)
   ..                   

      先行する計算での各反復での(Bi)CG法のパラメーター :math:`\beta`.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: z_seed ! ("CG_R_restart", "cg_c_restart"の場合)
      COMPLEX(8),INTENT(IN) :: z_seed ! (それ以外)
   ..                   

      先行する計算でのシードシフト.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: r_l_save(nl,iter_old) ! ("CG_R_restart"の場合)
      COMPLEX(8),INTENT(IN) :: r_l_save(nl,iter_old) ! (それ以外)
   ..                   

      先行する計算での各反復での射影された残差ベクトル.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: v4(ndim) ! ("CG_R_restart" の場合)
      COMPLEX(8),INTENT(IN) :: v4(ndim) ! (それ以外)
   ..
   
      ``BiCG_restart`` の場合のみ使用.
      先行する計算での最後の影の残差ベクトル.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: v14(ndim) ! ("CG_R_restart" の場合)
      COMPLEX(8),INTENT(IN) :: v14(ndim) ! (それ以外)
   ..

      ``BiCG_restart`` の場合のみ使用.
      先行する計算での最後から2番目の影の残差ベクトル.

.. _update:
      
\*_update
~~~~~~~~~

ループ内で行列ベクトル積と交互に呼ばれて解を更新する.

構文

   Fortran (シリアル/OpenMPI版)

   .. code-block:: fortran

      CALL komega_cg_r_update(v12, v2, x, r_l, status)
      CALL komega_cg_c_update(v12, v2, x, r_l, status)
      CALL komega_cocg_update(v12, v2, x, r_l, status)
      CALL komega_bicg_update(v12, v2, v14, v4, x, r_l, status)

   Fortran (MPI/ハイブリッド並列版)

   .. code-block:: fortran

      CALL pkomega_cg_r_update(v12, v2, x, r_l, status)
      CALL pkomega_cg_c_update(v12, v2, x, r_l, status)
      CALL pkomega_cocg_update(v12, v2, x, r_l, status)
      CALL pkomega_bicg_update(v12, v2, v14, v4, x, r_l, status)

   C/C++ (シリアル/OpenMPI版)

   .. code-block:: c

      komega_cg_r_update(v12, v2, x, r_l, status);
      komega_cg_c_update(v12, v2, x, r_l, status);
      komega_cocg_update(v12, v2, x, r_l, status);
      komega_bicg_update(v12, v2, v14, v4, x, r_l, status);

   C/C++ (MPI/ハイブリッド並列版)

   .. code-block:: c

      pkomega_cg_r_update(v12, v2, x, r_l, status);
      pkomega_cg_c_update(v12, v2, x, r_l, status);
      pkomega_cocg_update(v12, v2, x, r_l, status);
      pkomega_bicg_update(v12, v2, v14, v4, x, r_l, status);

   パラメーター

   .. code-block:: fortran

      REAL(8),INTENT(INOUT) :: v12(ndim) ! ("CG_R_update" の場合)
      COMPLEX(8),INTENT(INOUT) :: v12(ndim) ! (それ以外)
   ..

      入力は残差ベクトル( ``v2``)と行列の積. 出力は,
      更新された残差ベクトルの2-ノルムが,
      先頭の要素に格納される(これは収束の具合を表示して調べる時などに用いる).

   .. code-block:: fortran

      REAL(8),INTENT(INOUT) :: v2(ndim) ! ("CG_R_update" の場合)
      COMPLEX(8),INTENT(INOUT) :: v2(ndim) ! (それ以外)
   ..
   
      入力は残差ベクトル.
      出力は更新された残差ベクトル.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: v14(ndim) ! ("CG_R_update" の場合)
      COMPLEX(8),INTENT(IN) :: v14(ndim) ! (それ以外)
   ..

      影の残差ベクトル( ``v4``)と行列の積.

   .. code-block:: fortran

      REAL(8),INTENT(INOUT) :: v4(ndim) ! ("CG_R_update" の場合)
      COMPLEX(8),INTENT(INOUT) :: v4(ndim) ! (それ以外)
   ..

      入力は影の残差ベクトル.
      出力は更新された影の残差ベクトル.

   .. code-block:: fortran

      INTEGER,INTENT(OUT) :: status(3)
   ..
   
      エラーコードを返す.

      第一成分( ``status(1)``)
      
         解が収束した場合,
         もしくは計算が破綻した場合には現在の総反復回数に
         マイナスが付いた値が返される.
         それ以外の場合には現在の総反復回数(マイナスが付かない)が返される.
         ``status(1)`` が正の値の時のみ反復を続行できる.
         それ以外の場合は反復を進めても有意な結果は得られない.

      第二成分( ``status(2)``)
      
         :ref:`init` ルーチンで, ``itermax`` を有限にして,
         かつ ``itermax`` 回の反復で
         収束に達しなかった場合には ``1`` が返される.
         :math:`\alpha` が発散した場合には ``2`` が返される.
         :math:`\pi_{\rm seed}` が0にになった場合には ``3`` が返される.
         ``COCG_update`` もしくは ``BiCG_update`` で,
         残差ベクトルと影の残差ベクトルが直交した場合には ``4`` が返される.
         それ以外の場合には ``0`` が返される.

      第三成分( ``status(3)``)
      
         シード点のindexが返される.

.. _getcoef:
         
\*_getcoef
~~~~~~~~~~

後でリスタートをするときに必要な係数を取得する.
このルーチンを呼び出すためには,
:ref:`init` ルーチンで ``itermax`` を ``0`` 以外の値にしておく必要がある.
     
また, このルーチンで使われる総反復回数 (``iter_old``) は :ref:`update` の出力 ``status``
を用いて次のように計算される.

.. code-block:: fortran

   iter_old = ABS(status(1))

構文

   Fortran (シリアル/OpenMP版)

   .. code-block:: fortran

       CALL komega_cg_r_getcoef(alpha_save, beta_save, z_seed, r_l_save)
       CALL komega_cg_c_getcoef(alpha_save, beta_save, z_seed, r_l_save)
       CALL komega_cocg_getcoef(alpha_save, beta_save, z_seed, r_l_save)
       CALL komega_bicg_getcoef(alpha_save, beta_save, z_seed, r_l_save)

   Fortran (MPI/ハイブリッド並列版)

   .. code-block:: fortran

       CALL pkomega_cg_r_getcoef(alpha_save, beta_save, z_seed, r_l_save)
       CALL pkomega_cg_c_getcoef(alpha_save, beta_save, z_seed, r_l_save)
       CALL pkomega_cocg_getcoef(alpha_save, beta_save, z_seed, r_l_save)
       CALL pkomega_bicg_getcoef(alpha_save, beta_save, z_seed, r_l_save)

   C/C++ (シリアル/OpenMP版)

   .. code-block:: c

       komega_cg_r_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
       komega_cg_c_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
       komega_cocg_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
       komega_bicg_getcoef(alpha_save, beta_save, &z_seed, r_l_save);

   C/C++ (MPI/ハイブリッド並列版)

   .. code-block:: c

       pkomega_cg_r_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
       pkomega_cg_c_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
       pkomega_cocg_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
       pkomega_bicg_getcoef(alpha_save, beta_save, &z_seed, r_l_save);

パラメーター

   .. code-block:: fortran

      REAL(8),INTENT(OUT) :: alpha_save(iter_old) ! ("CG_R_restart", "cg_c_restart"の場合)
      COMPLEX(8),INTENT(OUT) :: alpha_save(iter_old) ! (それ以外)
   ..
   
      各反復での(Bi)CG法のパラメーター :math:`\alpha`.

   .. code-block:: fortran

      REAL(8),INTENT(OUT) :: beta_save(iter_old) ! ("CG_R_restart", "cg_c_restart"の場合)
      COMPLEX(8),INTENT(OUT) :: beta_save(iter_old) ! (それ以外)
   ..                   

      各反復での(Bi)CG法のパラメーター :math:`\beta`.

   .. code-block:: fortran

      REAL(8),INTENT(OUT) :: z_seed ! ("CG_R_restart", "cg_c_restart"の場合)
      COMPLEX(8),INTENT(OUT) :: z_seed ! (それ以外)
   ..                   

      シードシフト.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: r_l_save(nl,iter_old) ! ("CG_R_restart"の場合)
      COMPLEX(8),INTENT(IN) :: r_l_save(nl,iter_old) ! (それ以外)
   ..                   

      各反復での射影された残差ベクトル.

.. _getvec:
      
\*_getvec
~~~~~~~~~

後でリスタートをするときに必要な残差ベクトルを取得する.
このルーチンを呼び出すためには,
:ref:`init` ルーチンで ``itermax`` を ``0`` 以外の値にしておく必要がある.

構文

   Fortran (シリアル/OpenMP版)

   .. code-block:: fortran

       CALL komega_cg_r_getvec(r_old)
       CALL komega_cg_c_getvec(r_old)
       CALL komega_cocg_getvec(r_old)
       CALL komega_bicg_getvec(r_old, r_tilde_old)

   Fortran (MPI/ハイブリッド並列版)

   .. code-block:: fortran

       CALL pkomega_cg_r_getvec(r_old)
       CALL pkomega_cg_c_getvec(r_old)
       CALL pkomega_cocg_getvec(r_old)
       CALL pkomega_bicg_getvec(r_old, r_tilde_old)

   C/C++ (シリアル/OpenMP版)

   .. code-block:: c

       komega_cg_r_getvec(r_old);
       komega_cg_c_getvec(r_old);
       komega_cocg_getvec(r_old);
       komega_bicg_getvec(r_old, r_tilde_old);

   C/C++ (MPI/ハイブリッド並列版)

   .. code-block:: c

       pkomega_cg_r_getvec(r_old);
       pkomega_cg_c_getvec(r_old);
       pkomega_cocg_getvec(r_old);
       pkomega_bicg_getvec(r_old, r_tilde_old);

パラメーター

   .. code-block:: fortran

      REAL(8),INTENT(OUT) :: r_old(ndim) ! ("CG_R_getvec" の場合)
      COMPLEX(8),INTENT(OUT) :: r_old(ndim) ! (それ以外)
   ..

      先行する計算での最後から2番目の残差ベクトル.

   .. code-block:: fortran

      COMPLEX(8),INTENT(OUT) :: r_tilde_old(ndim)
   ..

      ``BiCG_getvec`` の場合のみ使用.
      先行する計算での最後から2番目の影の残差ベクトル.

\*_getresidual
~~~~~~~~~~~~~~

各シフト点での残差ベクトルの2-ノルムを取得する. 
このルーチンは :ref:`init` と :ref:`finalize` の間の
任意の場所で呼び出すことが出来る. また,
いつ何回呼び出しても最終的な計算結果には影響を与えない.

構文

   Fortran (シリアル/OpenMP版)

   .. code-block:: fortran

       CALL komega_cg_r_getresidual(res)
       CALL komega_cg_c_getresidual(res)
       CALL komega_cocg_getresidual(res)
       CALL komega_bicg_getresidual(res)

   Fortran (MPI/ハイブリッド並列版)

   .. code-block:: fortran

       CALL pkomega_cg_r_getresidual(res)
       CALL pkomega_cg_c_getresidual(res)
       CALL pkomega_cocg_getresidual(res)
       CALL pkomega_bicg_getresidual(res)

   C/C++ (シリアル/OpenMP版)

   .. code-block:: c

       komega_cg_r_getresidual(res);
       komega_cg_c_getresidual(res);
       komega_cocg_getresidual(res);
       komega_bicg_getresidual(res);

   C/C++ (MPI/ハイブリッド並列版)

   .. code-block:: c

       pkomega_cg_r_getresidual(res);
       pkomega_cg_c_getresidual(res);
       pkomega_cocg_getresidual(res);
       pkomega_bicg_getresidual(res);

パラメーター

   .. code-block:: fortran

      COMPLEX(8),INTENT(OUT) :: res(nz)
   ..

      各シフト点での残差ベクトルの2-ノルム.

.. _finalize:
      
\*_finalize
~~~~~~~~~~~

ライブラリ内部で割りつけた配列のメモリを解放する.

構文

   Fortran (シリアル/OpenMP版)

   .. code-block:: fortran

       CALL komega_cg_r_finalize()
       CALL komega_cg_c_finalize()
       CALL komega_cocg_finalize()
       CALL komega_bicg_finalize()

   Fortran (MPI/ハイブリッド並列版)

   .. code-block:: fortran

       CALL pkomega_cg_r_finalize()
       CALL pkomega_cg_c_finalize()
       CALL pkomega_cocg_finalize()
       CALL pkomega_bicg_finalize()

   C/C++ (シリアル/OpenMP版)

   .. code-block:: c

       komega_cg_r_finalize();
       komega_cg_c_finalize();
       komega_cocg_finalize();
       komega_bicg_finalize();

   C/C++ (MPI/ハイブリッド並列版)

   .. code-block:: c

       pkomega_cg_r_finalize();
       pkomega_cg_c_finalize();
       pkomega_cocg_finalize();
       pkomega_bicg_finalize();

Shifted BiCGライブラリを使用したソースコードの例
------------------------------------------------

以下, 代表的な例としてShifted BiCGライブラリの場合の使用方法を記載する.

.. code-block:: fortran

   PROGRAM my_prog
     !
     USE komega_bicg, ONLY : komega_bicg_init, komega_bicg_restart, &
     &                       komega_bicg_update, komega_bicg_getcoef, &
     &                       komega_bicg_getvec, komega_bicg_finalize
     USE solve_cc_routines, ONLY : input_size, input_restart, &
     &                             projection, &
     &                             hamiltonian_prod, generate_system, &
     &                             output_restart, output_result
     !
     IMPLICIT NONE
     !
     INTEGER,SAVE :: &
     & ndim,    & ! Size of Hilvert space
     & nz,      & ! Number of frequencies
     & nl,      & ! Number of Left vector
     & itermax, & ! Max. number of iteraction
     & iter_old   ! Number of iteraction of previous run
     !
     REAL(8),SAVE :: &
     & threshold ! Convergence Threshold
     !
     COMPLEX(8),SAVE :: &
     & z_seed ! Seed frequency
     !
     COMPLEX(8),ALLOCATABLE,SAVE :: &
     & z(:)         ! (nz): Frequency
     !
     COMPLEX(8),ALLOCATABLE,SAVE :: &
     & ham(:,:), &
     & rhs(:), &
     & v12(:), v2(:), & ! (ndim): Working vector
     & v14(:), v4(:), & ! (ndim): Working vector
     & r_l(:), & ! (nl) : Projeccted residual vector 
     & x(:,:) ! (nl,nz) : Projected result 
     !
     ! Variables for Restart
     !
     COMPLEX(8),ALLOCATABLE,SAVE :: &
     & alpha(:), beta(:) ! (iter_old) 
     !
     COMPLEX(8),ALLOCATABLE,SAVE :: &
     & r_l_save(:,:) ! (nl,iter_old) Projected residual vectors
     !
     ! Variables for Restart
     !
     INTEGER :: &
     & iter,    & ! Counter for Iteration
     & status(3)
     !
     LOGICAL :: &
     & restart_in, & ! If .TRUE., sestart from the previous result
     & restart_out   ! If .TRUE., save datas for the next run
     !
     ! Input Size of vectors, numerical conditions
     !
     CALL input_size(ndim,nl,nz)
     CALL input_condition(itermax,threshold,restart_in,restart_out)
     !
     ALLOCATE(v12(ndim), v2(ndim), v14(ndim), v4(ndim), r_l(nl), &
     &        x(nl,nz), z(nz), ham(ndim,ndim), rhs(ndim))
     !
     CALL generate_system(ndim, ham, rhs, z)
     !
     WRITE(*,*)
     WRITE(*,*) "#####  CG Initialization  #####"
     WRITE(*,*)
     !
     IF(restart_in) THEN
       !
       CALL input_restart(iter_old, zseed, alpha, beta, r_l_save)
       !
       IF(restart_out) THEN
          CALL komega_bicg_restart( &
          &    ndim, nl, nz, x, z, itermax, threshold, &
          &    status, iter_old, v2, v12, v4, v14, alpha, &
          &    beta, z_seed, r_l_save)
       ELSE
          CALL komega_bicg_restart( &
          &    ndim, nl, nz, x, z, 0, threshold, &
          &    status, iter_old, v2, v12, v4, v14, alpha, &
          &    beta, z_seed, r_l_save)
       END IF
       !
       ! These vectors were saved in BiCG routine
       !
       DEALLOCATE(alpha, beta, r_l_save)
       !
       IF(status(1) /= 0) GOTO 10
       !
     ELSE
        !
        ! Generate Right Hand Side Vector
        !
        v2(1:ndim) = rhs(1:ndim)
        v4(1:ndim) = CONJG(v2(1:ndim))
        !v4(1:ndim) = v2(1:ndim)
        !
        IF(restart_out) THEN
           CALL komega_bicg_init(ndim, nl, nz, x, z, termax, threshold)
        ELSE
           CALL komega_bicg_init(ndim, nl, nz, x, z, 0, threshold)
        END IF
        !
     END IF
     !
     ! BiCG Loop
     !
     WRITE(*,*)
     WRITE(*,*) "#####  CG Iteration  #####"
     WRITE(*,*)
     !
     DO iter = 1, itermax
        !
        ! Projection of Residual vector into the space
        ! spaned by left vectors
        !
        r_l(1:nl) = projection(v2(1:nl))
        !
        ! Matrix-vector product
        !
        CALL hamiltonian_prod(Ham, v2, v12)
        CALL hamiltonian_prod(Ham, v4, v14)
        !
        ! Update result x with BiCG
        !
        CALL komega_bicg_update(v12, v2, v14, v4, x, r_l, status)
        !
        WRITE(*,'(a,i,a,3i,a,e15.5)') "lopp : ", iter, &
        &                             ", status : ", status(1:3), &
        &                             ", Res. : ", DBLE(v12(1))
        IF(status(1) < 0) EXIT
        !
     END DO
     !
     IF(status(2) == 0) THEN
        WRITE(*,*) "  Converged in iteration ", ABS(status(1))
     ELSE IF(status(2) == 1) THEN
        WRITE(*,*) "  Not Converged in iteration ", ABS(status(1))
     ELSE IF(status(2) == 2) THEN
        WRITE(*,*) "  Alpha becomes infinity", ABS(status(1))
     ELSE IF(status(2) == 3) THEN
        WRITE(*,*) "  Pi_seed becomes zero", ABS(status(1))
     ELSE IF(status(2) == 4) THEN
     WRITE(*,*) "  Residual & Shadow residual are orthogonal", &
     &          ABS(status(1))
     END IF
     !
     ! Total number of iteration
     !
     iter_old = ABS(status(1))
     !
     ! Get these vectors for restart in the Next run
     !
     IF(restart_out) THEN
        !
        ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl,iter_old))
        !
        CALL komega_bicg_getcoef(alpha, beta, z_seed, r_l_save)
        CALL komega_bicg_getvec(v12,v14)
        !
        CALL output_restart(iter_old, z_seed, alpha, beta, &
        &                   r_l_save, v12, v14)
        !
        DEALLOCATE(alpha, beta, r_l_save)
        !     
     END IF
     !
   10 CONTINUE
     !
     ! Deallocate all intrinsic vectors
     !
     CALL komega_bicg_finalize()
     !
     ! Output to a file
     !
     CALL output_result(nl, nz, z, x, r_l)
     !
     DEALLOCATE(v12, v2, v14, v4, r_l, x, z)
     !
     WRITE(*,*)
     WRITE(*,*) "#####  Done  #####"
     WRITE(*,*)
     !
   END PROGRAM my_prog

