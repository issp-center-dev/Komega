プログラム内でのライブラリの動作イメージ
========================================

以下では :math:`N_R` のループは省略する(各右辺ベクトルごとに同じ事をすればいいので).
また :math:`G_{i j}(z_k)` の代わりに :math:`N_z` 個の :math:`N_L` 次元の解ベクトル :math:`{\bf x}_{k}` を求める.

注意事項.

-  ``komega_????_init`` を呼び出す前にサイズ :math:`N_H` のベクトルを2本(BiCGの時には4本)Allocateしておく.

-  ハミルトニアン-ベクトル積を行う部分はあらかじめ作成しておく.

-  解ベクトルをAllocateしておく. ただし,
   これは必ずしも :math:`N_H \times N_z` である必要はない.
   実際前節の場合は :math:`N_L \times N_z` である.この時(双)共役勾配ベクトル :math:`{{\bf p}_k}` も
   :math:`N_z` 個の :math:`N_L` 次元のベクトルである.
   :math:`N_H` 次元の残差ベクトルを :math:`N_L` 次元へ変換するところはユーザーが行う.

   .. math::

      \begin{aligned}
          {\bf r}^{\rm L} = {\hat P}^\dagger {\boldsymbol r}, \qquad
          {\hat P} \equiv ({\boldsymbol \varphi}_1, \cdots, {\boldsymbol \varphi}_{N_L})
        \end{aligned}

-  ``komega_????_update`` の出力 ``status`` の第一成分が負の値になった場合には,
   解が収束した, 若しくは破たんしたことを表す.
   したがって ``status(1) < 0`` でループを抜けるようにしておく.

-  ``komega_????_update`` 内での収束判定には,
   シード点での残差ベクトルの2-ノルムが使われる. すなわち,
   すべてのシフト点での残差ベクトルの2-ノルムが ``threshold`` を下回った時に収束したと見做される.

-  各iterationでの :math:`\alpha, \beta, {\bf r}^{\rm L}` を保存しておき,
   あとで利用する場合には最大反復回数``itermax`` を設定する.

ライブラリの各ルーチンの名前は次の通りである.

-  ``komega_BiCG_init``, ``komega_COCG_init``, ``komega_CG_C_init``,
   ``komega_CG_R_init``

   ライブラリ内部で使う(ユーザーの目に触れない)変数のAllocateや初期値設定を行う.

-  ``komega_BiCG_update``, ``komega_COCG_update``,
   ``komega_CG_C_update``, ``komega_CG_R_update``

   Iteration の中で呼び出される. 解ベクトル群の更新等を行う.

-  ``komega_BiCG_finalize``, ``komega_COCG_finalize``,
   ``komega_CG_C_finalize``, ``komega_CG_R_finalize``

   Allocateしたライブラリ内部ベクトルを開放する.

-  ``komega_BiCG_getcoef``, ``komega_COCG_getcoef``,
   ``komega_CG_C_getcoef``, ``komega_CG_R_getcoef``

   各iterationで保存しておいた :math:`\alpha`, :math:`\beta`,
   :math:`z_{\rm seed}`, :math:`{\bf r}^{\rm L}` を取り出す.

-  ``komega_BiCG_getvec``, ``komega_COCG_getvec``,
   ``komega_CG_C_getvec``, ``komega_CG_R_getvec``

   :math:`{\boldsymbol r}`, :math:`{\boldsymbol r}^{\rm old}`,
   :math:`{\tilde {\boldsymbol r}}`,
   :math:`{\tilde {\boldsymbol r}}^{\rm old}` を取り出す.

-  ``komega_BiCG_restart``, ``komega_COCG_restart``,
   ``komega_CG_C_restart``, ``komega_CG_R_restart``

   保存しておいた :math:`\alpha` 等を用いて,
   新規の :math:`z` での計算を行う.
   :math:`{\boldsymbol r}` 等も有る場合には ``komega_BiCG_init``,
   ``komega_COCG_init``, ``komega_CG_C_init``, ``komega_CG_R_init``
   の代わりに用いてリスタートすることもできる.

Shifted BiCGライブラリの動作イメージ
------------------------------------

Allocate :math:`{\boldsymbol v}_{1 2}`, :math:`{\boldsymbol v}_{1 3}`,
:math:`{\boldsymbol v}_2`, :math:`{\boldsymbol v}_3`,
:math:`\{{\bf x}_k\}, {\bf r}^{\rm L}`
:math:`{\boldsymbol v}_2 = {\boldsymbol \varphi_j}`

``komega_BiCG_init(N_H, N_L, N_z, x, z, itermax, threshold)`` start

   Allocate :math:`{\boldsymbol v}_3`, :math:`{\boldsymbol v}_5`,
   :math:`\{\pi_k\}`, :math:`\{\pi_k^{\rm old}\}`, :math:`\{{\bf p}_k\}`

   Copy :math:`\{z_k\}`

   ``itermax`` :math:`\neq` ``0`` ならば :math:`\alpha`,
   :math:`\beta`, :math:`{\bf r}^{\rm L}` を保存する配列を確保する.

   :math:`{\boldsymbol v}_4 = {\boldsymbol v}_2^*` (任意),
   :math:`{\boldsymbol v}_3 = {\boldsymbol v}_5 = {\bf 0}`,

   :math:`{\bf p}_{k} = {\bf x}_k = {\bf 0}(k=1 \cdots N_z),\; \pi_k=\pi_k^{\rm old} = 1(k=1 \cdots N_z)`

   :math:`\rho = \infty,\; \alpha = 1,\; z_{\rm seed}=0`

   ( :math:`{\boldsymbol v}_2 \equiv {\boldsymbol r}`,
   :math:`{\boldsymbol v}_3 \equiv {\boldsymbol r}^{\rm old}`,
   :math:`{\boldsymbol v}_4 \equiv {\tilde {\boldsymbol r}}`,
   :math:`{\boldsymbol v}_5 \equiv {\tilde {\boldsymbol r}}^{\rm old}`. )

``komega_BiCG_init`` finish

do iteration

   :math:`{\bf r}^{\rm L} = {\hat P}^\dagger {\boldsymbol v}_2`

   :math:`{\boldsymbol v}_{1 2} = {\hat H} {\boldsymbol v}_2`,
   :math:`{\boldsymbol v}_{1 4} = {\hat H} {\boldsymbol v}_4`
   [ :math:`({\boldsymbol v}_{1 2}, {\boldsymbol v}_{1 4}) = {\hat H} ({\boldsymbol v}_2, {\boldsymbol v}_4)` とも書ける]

   ``komega_BiCG_update(v_12, v_2, v_14, v_4, x, r_small, status)`` start

      :math:`\circ` シード方程式

      :math:`\rho^{\rm old} = \rho,\; \rho = {\boldsymbol v}_4^* \cdot {\boldsymbol v}_2`

      :math:`\beta = \rho / \rho^{\rm old}`

      :math:`{\boldsymbol v}_{1 2} = z_{\rm seed} {\boldsymbol v}_2 - {\boldsymbol v}_{1 2}`,
      :math:`{\boldsymbol v}_{1 4} = z_{\rm seed}^* {\boldsymbol v}_4 - {\boldsymbol v}_{1 4}`

      :math:`\alpha^{\rm old} = \alpha,\; \alpha = \frac{\rho}{{\boldsymbol v}_3^* \cdot {\boldsymbol v}_{1 2} - \beta \rho / \alpha }`

      :math:`\circ` シフト方程式

      do :math:`k = 1 \cdots N_z`

         :math:`\pi_k^{\rm new} = [1+\alpha(z_k-z_{\rm seed})]\pi_k - \frac{\alpha \beta}{\alpha^{\rm old}}(\pi_k^{\rm old} - \pi_k)`

         :math:`{\bf p}_{k} = \frac{1}{\pi_k} {\bf r}^{\rm L} + \frac{\pi^{\rm old}_k \pi^{\rm old}_k}{\pi_k \pi_k} \beta {\bf p}_{k}`

         :math:`{\bf x}_{k} = {\bf x}_{k} + \frac{\pi_k}{\pi_k^{\rm new}} \alpha {\bf p}_{k}`

         :math:`\pi_k^{\rm old} = \pi_k`, :math:`\pi_k = \pi_k^{\rm new}`

      end do :math:`k`

      :math:`{\boldsymbol v}_{1 2} = \left( 1 + \frac{\alpha \beta}{\alpha^{\rm old}} \right) {\boldsymbol v}_2 - \alpha {\boldsymbol v}_{1 2} - \frac{\alpha \beta}{\alpha^{\rm old}} {\boldsymbol v}_3`,
      :math:`{\boldsymbol v}_3 = {\boldsymbol v}_2,\; {\boldsymbol v}_2 = {\boldsymbol v}_{1 2}`

      :math:`{\boldsymbol v}_{1 4} = \left( 1 + \frac{\alpha^* \beta^*}{\alpha^{{\rm old}*}} \right) {\boldsymbol v}_4 - \alpha^* {\boldsymbol v}_{1 4} - \frac{\alpha^* \beta^*}{\alpha^{{\rm old} *}} {\boldsymbol v}_5`,
      :math:`{\boldsymbol v}_5 = {\boldsymbol v}_4,\; {\boldsymbol v}_4 = {\boldsymbol v}_{1 4}`

      :math:`\circ` Seed switch

      :math:`|\pi_k|` が最も小さい :math:`k` を探す.
      :math:`\rightarrow z_{\rm seed},\; \pi_{\rm seed},\; \pi_{\rm seed}^{\rm old}`

      :math:`{\boldsymbol v}_2 = {\boldsymbol v}_2 / \pi_{\rm seed}`,
      :math:`{\boldsymbol v}_3 = {\boldsymbol v}_3 / \pi_{\rm seed}^{\rm old}`,
      :math:`{\boldsymbol v}_4 = {\boldsymbol v}_4 / \pi_{\rm seed}^{*}`,
      :math:`{\boldsymbol v}_5 = {\boldsymbol v}_5 / \pi_{\rm seed}^{\rm old *}`

      :math:`\alpha = (\pi_{\rm seed}^{\rm old} / \pi_{\rm seed}) \alpha`,
      :math:`\rho = \rho / (\pi_{\rm seed}^{\rm old} \pi_{\rm seed}^{\rm old})`

      :math:`\{\pi_k = \pi_k / \pi_{\rm seed},\; \pi_k^{\rm old} = \pi_k^{\rm old} / \pi_{\rm seed}^{\rm old}\}`

   ``komega_BiCG_update`` finish

   if(status(1) < 0 (これは:math:`|{\boldsymbol v}_2| <` Threshold となった事を意味する)) exit

end do iteration

``komega_BiCG_finalize`` start

   Deallocate :math:`{\boldsymbol v}_4`, :math:`{\boldsymbol v}_5`,
   :math:`\{\pi_k\}`, :math:`\{\pi_k^{\rm old}\}`, :math:`\{{\bf p}_k\}`

``komega_BiCG_finalize`` finish

Shifted COCGライブラリの動作イメージ
------------------------------------

Allocate :math:`{\boldsymbol v}_1`, :math:`{\boldsymbol v}_2`,
:math:`\{{\bf x}_k\}, {\bf r}^{\rm L}`
:math:`{\boldsymbol v}_2 = {\boldsymbol \varphi_j}`

``COCG_init(N_H, N_L, N_z, x, z, itermax, threshold)`` start

   Allocate :math:`{\boldsymbol v}_3`, :math:`\{\pi_k\}`,
   :math:`\{\pi_k^{\rm old}\}`, :math:`\{{\bf p}_k\}`

   Copy :math:`\{z_k\}`

   ``itermax`` :math:`\neq` ``0`` ならば :math:`\alpha`,
   :math:`\beta`, :math:`{\bf r}^{\rm L}` を保存する配列を確保する.

   :math:`{\boldsymbol v}_3 = {\bf 0}`,

   :math:`{\bf p}_{k} = {\bf x}_k = {\bf 0}(k=1 \cdots N_z),\; \pi_k=\pi_k^{\rm old} = 1(k=1 \cdots N_z)`

   :math:`\rho = \infty,\; \alpha = 1,\; \beta=0,\; z_{\rm seed}=0`

   ( :math:`{\boldsymbol v}_2 \equiv {\boldsymbol r}`,
   :math:`{\boldsymbol v}_3 \equiv {\boldsymbol r}^{\rm old}`. )
         
``COCG_init`` finish

do iteration

   :math:`{\bf r}^{\rm L} = {\hat P}^\dagger {\boldsymbol v}_2`

   :math:`{\boldsymbol v}_1 = {\hat H} {\boldsymbol v}_2`

   ``COCG_update(v_1, v_2, x, r_small, status)`` start

      :math:`\circ` シード方程式

      :math:`\rho^{\rm old} = \rho,\; \rho = {\boldsymbol v}_2 \cdot {\boldsymbol v}_2`

      :math:`\beta = \rho / \rho^{\rm old}`

      :math:`{\boldsymbol v}_1 = z_{\rm seed} {\boldsymbol v}_2 - {\boldsymbol v}_1`

      :math:`\alpha^{\rm old} = \alpha,\; \alpha = \frac{\rho}{{\boldsymbol v}_2 \cdot {\boldsymbol v}_1 - \beta \rho / \alpha }`

      :math:`\circ` シフト方程式

      do :math:`k = 1 \cdots N_z`

         :math:`\pi_k^{\rm new} = [1+\alpha(z_k-z_{\rm seed})]\pi_k - \frac{\alpha \beta}{\alpha^{\rm old}}(\pi_k^{\rm old} - \pi_k)`

         :math:`{\bf p}_{k} = \frac{1}{\pi_k} {\bf r}^{\rm L} + \frac{\pi^{\rm old}_k \pi^{\rm old}_k}{\pi_k \pi_k} \beta {\bf p}_{k}`

         :math:`{\bf x}_{k} = {\bf x}_{k} + \frac{\pi_k}{\pi_k^{\rm new}} \alpha {\bf p}_{k}`

         :math:`\pi_k^{\rm old} = \pi_k,\; \pi_k = \pi_k^{\rm new}`

      end do :math:`k`

      :math:`{\boldsymbol v}_1 = \left( 1 + \frac{\alpha \beta}{\alpha^{\rm old}} \right) {\boldsymbol v}_2 - \alpha {\boldsymbol v}_1 - \frac{\alpha \beta}{\alpha^{\rm old}} {\boldsymbol v}_3`

      :math:`{\boldsymbol v}_3 = {\boldsymbol v}_2`,
      :math:`{\boldsymbol v}_2 = {\boldsymbol v}_1`

      :math:`\circ` Seed switch

      :math:`|\pi_k|` が最も小さい :math:`k` を探す.
      :math:`\rightarrow z_{\rm seed},\; \pi_{\rm seed},\; \pi_{\rm seed}^{\rm old}`

      :math:`{\boldsymbol v}_2 = {\boldsymbol v}_2 / \pi_{\rm seed}`,
      :math:`{\boldsymbol v}_3 = {\boldsymbol v}_3 / \pi_{\rm seed}^{\rm old}`

      :math:`\alpha = (\pi_{\rm seed}^{\rm old} / \pi_{\rm seed}) \alpha`,
      :math:`\rho = \rho / (\pi_{\rm seed}^{\rm old} \pi_{\rm seed}^{\rm old})`

      :math:`\{\pi_k = \pi_k / \pi_{\rm seed},\; \pi_k^{\rm old} = \pi_k^{\rm old} / \pi_{\rm seed}^{\rm old}\}`

   ``COCG_update`` finish

   if(status(1) < 0 (これは:math:`|{\boldsymbol v}_2| <` Threshold となった事を意味する)) exit

end do iteration

``COCG_finalize`` start

   Deallocate :math:`{\boldsymbol v}_3`, :math:`\{\pi_k\}`,
   :math:`\{\pi_k^{\rm old}\}`, :math:`\{{\bf p}_k\}`

``COCG_finalize`` finish

Shifted CGライブラリの動作イメージ
----------------------------------

COCGと同様.

