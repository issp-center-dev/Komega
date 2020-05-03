概要
====

本資料はISSP Math
Libraryの内の, Krylov部分空間法に基づくシフト線形方程式群ライブラリ
:math:`K\omega` を用いた Green関数計算用ミニアプリのマニュアルです. 
ライブラリに関する使用方法については, " :math:`K\omega` マニュアル"に記載しました. 

ソフトウェア概要
----------------

本ソフトウェアでは, Green関数

.. math::

   \begin{align}
   G_{i}(z) =
   \langle i | (z-{\hat H})^{-1}| i \rangle
   \equiv 
   {\boldsymbol \varphi}_i^{*} \cdot (z-{\hat H})^{-1} {\boldsymbol \varphi}_i
   \end{align}

の計算を行います. 
ここで :math:`| i \rangle` はベクトル, :math:`{\cal H}` はハミルトニアン, 
:math:`z` は複素数シフトを表します. 

なお :math:`{\cal H}` については, 

-  :math:`{\cal H}` をMatrixMarket形式の入力ファイルとして与えるモード

-  Heisenberg模型の :math:`{\cal H}` を内部で与えるモード

を用意します. 
またグリーン関数の計算では, :math:`{\hat H}` と :math:`z` が複素数もしくは実数かに応じ, 

-  :math:`{\hat H}` も :math:`z` も両方複素数の場合 : Shifted
   Bi-Conjugate Gradient(BiCG)法

-  :math:`{\hat H}` が実数で :math:`z` が複素数の場合 : Shifted
   Conjugate Orthogonal Conjugate Gradient(COCG)法

の手法を用意しています.

.. デフォルトでは入力条件および入力ファイルからどのタイプかを自動判定し計算が行われます. 

.. ただし, 各手法を手動で選択させ計算させることも可能となっています. 

