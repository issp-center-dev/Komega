アルゴリズム
============

このライブラリは,
:math:`{\hat H}` および :math:`z` が複素数であるか実数であるかに応じて,
次の4種類の計算をサポートする( :math:`{\hat H}` は複素数の場合はエルミート行列,
実数の場合は実対称行列).

-  :math:`{\hat H}` も :math:`z` も両方複素数の場合 : Shifted
   Bi-Conjugate Gradient(BiCG)法 [Flommer2003]_

-  :math:`{\hat H}` が実数で :math:`z` が複素数の場合 : Shifted
   Conjugate Orthogonal Conjugate Gradient(COCG)法 [Yamamoto2008]_

-  :math:`{\hat H}` が複素数で :math:`z` が実数の場合 : Shifted
   Conjugate Gradient(CG)法 (複素ベクトル)

-  :math:`{\hat H}` も :math:`z` も両方実数の場合 : Shifted Conjugate
   Gradient(CG)法 (実ベクトル)

いずれの場合も Seed switching [Yamamoto2008]_ を行う. 左ベクトルが :math:`N_L` 個,
右ベクトルが :math:`N_R` 個(典型的には1個)あるとする. 以下,
各手法のアルゴリズムを記載する.

.. [Flommer2003] A. Frommer, Computing **70**, 87 (2003).

.. [Yamamoto2008] S. Yamamoto, T. Sogabe, T. Hoshi, S.-L. Zhang, and T. Fujiwara, J. Phys. Soc. Jpn. **77**, 114713 (2008).

Seed switch 付き Shifted BiCG法
-------------------------------

:math:`G_{i j}(z_k) = 0 (i=1 \cdots N_L,\; j = 1 \cdots N_R,\; k=1 \cdots N_z)`

do :math:`j = 1 \cdots N_R`

   :math:`{\boldsymbol r} = {\boldsymbol \varphi_j}`,

   :math:`{\tilde {\boldsymbol r}} =` 任意,
   :math:`{\boldsymbol r}^{\rm old} = {\tilde {\boldsymbol r}}^{\rm old} = {\bf 0}`

   :math:`p_{i k} = 0(i=1 \cdots N_L,\; k=1 \cdots N_z),\; \pi_k=\pi_k^{\rm old} = 1(k=1 \cdots N_z)`

   :math:`\rho = \infty,\; \alpha = 1,\; z_{\rm seed}=0`

   do iteration

      :math:`\circ` シード方程式

      :math:`\rho^{\rm old} = \rho,\; \rho = {\tilde {\boldsymbol r}}^* \cdot {\boldsymbol r}`

      :math:`\beta = \rho / \rho^{\rm old}`

      :math:`{\boldsymbol q} = (z_{\rm seed} {\hat I} - {\hat H}){\boldsymbol r}`

      :math:`\alpha^{\rm old} = \alpha,\; \alpha = \frac{\rho}{{\tilde {\boldsymbol r}}^*\cdot{\boldsymbol q} - \beta \rho / \alpha }`

      :math:`\circ` シフト方程式

      do :math:`k = 1 \cdots N_z`

         :math:`\pi_k^{\rm new} = [1+\alpha(z_k-z_{\rm seed})]\pi_k - \frac{\alpha \beta}{\alpha^{\rm old}}(\pi_k^{\rm old} - \pi_k)`

         do :math:`i = 1 \cdots N_L`

            :math:`p_{i k} = \frac{1}{\pi_k} {\boldsymbol \varphi}_i^* \cdot {\boldsymbol r} + \frac{\pi^{\rm old}_k \pi^{\rm old}_k}{\pi_k \pi_k} \beta p_{i k}`

            :math:`G_{i j}(z_k) = G_{i j}(z_k) + \frac{\pi_k}{\pi_k^{\rm new}} \alpha p_{i k}`

            :math:`\pi_k^{\rm old} = \pi_k`, :math:`\pi_k = \pi_k^{\rm new}`

         end do :math:`i`

      end do :math:`k`

      :math:`{\boldsymbol q} = \left( 1 + \frac{\alpha \beta}{\alpha^{\rm old}} \right) {\boldsymbol r} - \alpha {\boldsymbol q} - \frac{\alpha \beta}{\alpha^{\rm old}} {\boldsymbol r}^{\rm old},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r},\; {\boldsymbol r} = {\boldsymbol q}`

      :math:`{\boldsymbol q} = (z_{\rm seed}^* {\hat I} - {\hat H}) {\tilde {\boldsymbol r}},\; {\boldsymbol q} = \left( 1 + \frac{\alpha^* \beta^*}{\alpha^{{\rm old}*}} \right) {\tilde {\boldsymbol r}} - \alpha^* {\boldsymbol q} - \frac{\alpha^* \beta^*}{\alpha^{{\rm old} *}} {\tilde {\boldsymbol r}}^{\rm old},\; {\tilde {\boldsymbol r}}^{\rm old} = {\tilde {\boldsymbol r}},\; {\tilde {\boldsymbol r}} = {\boldsymbol q}`

      :math:`\circ` Seed switch

      :math:`|\pi_k|` が最も小さい :math:`k` を探す. :math:`\rightarrow z_{\rm seed},\; \pi_{\rm seed},\; \pi_{\rm seed}^{\rm old}`

      :math:`{\boldsymbol r} = {\boldsymbol r} / \pi_{\rm seed},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r}^{\rm old} / \pi_{\rm seed}^{\rm old},\; {\tilde {\boldsymbol r}} = {\tilde {\boldsymbol r}} / \pi_{\rm seed}^*,\; {\tilde {\boldsymbol r}}^{\rm old} = {\tilde {\boldsymbol r}}^{\rm old} / \pi_{\rm seed}^{{\rm old}*}`

      :math:`\alpha = (\pi_{\rm seed}^{\rm old} / \pi_{\rm seed}) \alpha`, :math:`\rho = \rho / (\pi_{\rm seed}^{\rm old} \pi_{\rm seed}^{\rm old})`

      :math:`\{\pi_k = \pi_k / \pi_{\rm seed},\; \pi_k^{\rm old} = \pi_k^{\rm old} / \pi_{\rm seed}^{\rm old}\}`

      if( :math:`|{\boldsymbol r}| <` Threshold) exit

   end do iteration

end do :math:`j`

Seed switch 付き Shifted COCG法
-------------------------------

BiCGのアルゴリズムで,
:math:`{\tilde {\boldsymbol r}} = {\boldsymbol r}^*,\; {\tilde {\boldsymbol r}}^{\rm old} = {\boldsymbol r}^{{\rm old}*}` とすると得られる.

:math:`G_{i j}(z_k) = 0 (i=1 \cdots N_L,\; j = 1 \cdots N_R,\; k=1 \cdots N_z)`

do :math:`j = 1 \cdots N_R`

   :math:`{\boldsymbol r} = {\boldsymbol \varphi_j}`, :math:`{\boldsymbol r}^{\rm old} = {\bf 0}`

   :math:`p_{i k} = 0(i=1 \cdots N_L,\; k=1 \cdots N_z),\; \pi_k=\pi_k^{\rm old} = 1(k=1 \cdots N_z)`

   :math:`\rho = \infty,\; \alpha = 1,\; z_{\rm seed}=0`

   do iteration

      :math:`\circ` シード方程式

      :math:`\rho^{\rm old} = \rho,\; \rho = {\boldsymbol r} \cdot {\boldsymbol r}`

      :math:`\beta = \rho / \rho^{\rm old}`

      :math:`{\boldsymbol q} = (z_{\rm seed} {\hat I} - {\hat H}){\boldsymbol r}`

      :math:`\alpha^{\rm old} = \alpha,\; \alpha = \frac{\rho}{{\boldsymbol r}\cdot{\boldsymbol q} - \beta \rho / \alpha }`

      :math:`\circ` シフト方程式

      do :math:`k = 1 \cdots N_z`

         :math:`\pi_k^{\rm new} = [1+\alpha(z_k-z_{\rm seed})]\pi_k - \frac{\alpha \beta}{\alpha^{\rm old}}(\pi_k^{\rm old} - \pi_k)`

         do :math:`i = 1 \cdots N_L`

            :math:`p_{i k} = \frac{1}{\pi_k} {\boldsymbol \varphi}_i^* \cdot {\boldsymbol r} + \frac{\pi^{\rm old}_k \pi^{\rm old}_k}{\pi_k \pi_k} \beta p_{i k}`

            :math:`G_{i j}(z_k) = G_{i j}(z_k) + \frac{\pi_k}{\pi_k^{\rm new}} \alpha p_{i k}`

            :math:`\pi_k^{\rm old} = \pi_k`, :math:`\pi_k = \pi_k^{\rm new}`

         end do :math:`i`

      end do :math:`k`

      :math:`{\boldsymbol q} = \left( 1 + \frac{\alpha \beta}{\alpha^{\rm old}} \right) {\boldsymbol r} - \alpha {\boldsymbol q} - \frac{\alpha \beta}{\alpha^{\rm old}} {\boldsymbol r}^{\rm old},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r},\; {\boldsymbol r} = {\boldsymbol q}`

      :math:`\circ` Seed switch

      :math:`|\pi_k|` が最も小さい :math:`k` を探す. :math:`\rightarrow z_{\rm seed},\; \pi_{\rm seed},\; \pi_{\rm seed}^{\rm old}`
                  
      :math:`{\boldsymbol r} = {\boldsymbol r} / \pi_{\rm seed},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r}^{\rm old} / \pi_{\rm seed}^{\rm old}`

      :math:`\alpha = (\pi_{\rm seed}^{\rm old} / \pi_{\rm seed}) \alpha`, :math:`\rho = \rho / (\pi_{\rm seed}^{\rm old} \pi_{\rm seed}^{\rm old})`

      :math:`\{\pi_k = \pi_k/\pi_{\rm seed},\; \pi_k^{\rm old} = \pi_k^{\rm old} / \pi_{\rm seed}^{\rm old}\}`

      if( :math:`|{\boldsymbol r}| <` Threshold) exit

   end do iteration

end do :math:`j`

Seed switch 付き Shifted CG法
-----------------------------

BiCGのアルゴリズムで,
:math:`{\tilde {\boldsymbol r}} = {\boldsymbol r},\; {\tilde {\boldsymbol r}}^{\rm old} = {\boldsymbol r}^{\rm old}` とすると得られる.

:math:`G_{i j}(z_k) = 0 (i=1 \cdots N_L,\; j = 1 \cdots N_R,\; k=1 \cdots N_z)`

do :math:`j = 1 \cdots N_R`

   :math:`{\boldsymbol r} = {\boldsymbol \varphi_j}`, :math:`{\boldsymbol r}^{\rm old} = {\bf 0}`

   :math:`p_{i k} = 0(i=1 \cdots N_L,\; k=1 \cdots N_z),\; \pi_k=\pi_k^{\rm old} = 1(k=1 \cdots N_z)`

   :math:`\rho = \infty,\; \alpha = 1,\; z_{\rm seed}=0`

   do iteration

      :math:`\circ` シード方程式

      :math:`\rho^{\rm old} = \rho,\; \rho = {\boldsymbol r}^* \cdot {\boldsymbol r}`

      :math:`\beta = \rho / \rho^{\rm old}`

      :math:`{\boldsymbol q} = (z_{\rm seed} {\hat I} - {\hat H}){\boldsymbol r}`

      :math:`\alpha^{\rm old} = \alpha,\; \alpha = \frac{\rho}{{\boldsymbol r}^* \cdot {\boldsymbol q} - \beta \rho / \alpha }`

      :math:`\circ` シフト方程式

      do :math:`k = 1 \cdots N_z`

         :math:`\pi_k^{\rm new} = [1+\alpha(z_k-z_{\rm seed})]\pi_k - \frac{\alpha \beta}{\alpha^{\rm old}}(\pi_k^{\rm old} - \pi_k)`

         do :math:`i = 1 \cdots N_L`

            :math:`p_{i k} = \frac{1}{\pi_k} {\boldsymbol \varphi}_i^* \cdot {\boldsymbol r} + \left(\frac{\pi^{\rm old}_k}{\pi_k } \right)^2 \beta p_{i k}`

            :math:`G_{i j}(z_k) = G_{i j}(z_k) + \frac{\pi_k}{\pi_k^{\rm new}} \alpha p_{i k}`

            :math:`\pi_k^{\rm old} = \pi_k`, :math:`\pi_k = \pi_k^{\rm new}`

         end do :math:`i`

      end do :math:`k`

      :math:`{\boldsymbol q} = \left( 1 + \frac{\alpha \beta}{\alpha^{\rm old}} \right) {\boldsymbol r} - \alpha {\boldsymbol q} - \frac{\alpha \beta}{\alpha^{\rm old}} {\boldsymbol r}^{\rm old},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r},\; {\boldsymbol r} = {\boldsymbol q}`

      :math:`\circ` Seed switch

      :math:`|\pi_k|` が最も小さい :math:`k` を探す. :math:`\rightarrow z_{\rm seed},\; \pi_{\rm seed},\; \pi_{\rm seed}^{\rm old}`

      :math:`{\boldsymbol r} = {\boldsymbol r} / \pi_{\rm seed},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r}^{\rm old} / \pi_{\rm seed}^{\rm old}`

      :math:`\alpha = (\pi_{\rm seed}^{\rm old} / \pi_{\rm seed}) \alpha`, :math:`\rho = \rho / {\pi_{\rm seed}^{\rm old}}^2`

      :math:`\{\pi_k = \pi_k/\pi_{\rm seed},\; \pi_k^{\rm old} = \pi_k^{\rm old}/\pi_{\rm seed}^{\rm old}\}`

      if( :math:`|{\boldsymbol r}| <` Threshold) exit

   end do iteration

end do :math:`j`

