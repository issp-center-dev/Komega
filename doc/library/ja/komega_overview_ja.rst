概要
====

本資料はISSP Math
Libraryの内の、Krylov部分空間法に基づくシフト線形方程式群ソルバーライブラリ
:math:`K\omega`\ に関するマニュアルである. 本ライブラリは,
(射影付き)シフト線形問題

.. math::

   \begin{align}
     G_{i j}(z) = \langle i | (z {\hat I} -{\hat H})^{-1}| j \rangle \equiv 
     {\boldsymbol \varphi}_i^{*} \cdot (z{\hat I}-{\hat H})^{-1} {\boldsymbol \varphi}_j
     \end{align}

を, Krylov部分空間法を用いて解くためのルーチンを提供する.
言語はfortranを用いる. また, BLASレベル1ルーチンを使用する.
