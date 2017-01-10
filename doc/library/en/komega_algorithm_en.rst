Algorithm
=========

This library provides the four kinds of numerical solvers.
The kind of solvers is selected under the condition whether the Hamiltonian
:math:`{\hat H}` and/or the frequency :math:`z` are complex or real number.
It is noted that :math:`{\hat H}` must be Hermitian (symmetric)
for complex (real) number.

-  (:math:`{\hat H}`, :math:`z` ) = (complex, complex):
   Shifted Bi-Conjugate Gradient(BiCG) method [Flommer2003]_

-  (:math:`{\hat H}`, :math:`z` ) = (real, complex):
   Shifted Conjugate Orthogonal Conjugate Gradient(COCG) method [Yamamoto2008]_

-  (:math:`{\hat H}`, :math:`z` ) = (complex, real):
   Shifted Conjugate Gradient(CG) method (using complex vector)

-  (:math:`{\hat H}`, :math:`z` ) = (real, real):
   Shifted Conjugate Gradient(CG) method (using real vector)

For above methods, seed switching [Yamamoto2008]_ is adopted.
Hereafter, the number of the left (right) side vector is
written as :math:`N_L` (:math:`N_R`).
The details of each algorithm are written as follows.

.. [Flommer2003] A. Frommer, Computing **70**, 87 (2003).

.. [Yamamoto2008] S. Yamamoto, T. Sogabe, T. Hoshi, S.-L. Zhang, and T. Fujiwara, J. Phys. Soc. Jpn. **77**, 114713 (2008).

Shifted BiCG method with seed switching technique
-------------------------------------------------

:math:`G_{i j}(z_k) = 0 (i=1 \cdots N_L,\; j = 1 \cdots N_R,\; k=1 \cdots N_z)`

do :math:`j = 1 \cdots N_R`

   :math:`{\boldsymbol r} = {\boldsymbol \varphi_j}`,

   :math:`{\tilde {\boldsymbol r}} =` an arbitrary vector,
   :math:`{\boldsymbol r}^{\rm old} = {\tilde {\boldsymbol r}}^{\rm old} = {\bf 0}`

   :math:`p_{i k} = 0(i=1 \cdots N_L,\; k=1 \cdots N_z),\; \pi_k=\pi_k^{\rm old} = 1(k=1 \cdots N_z)`

   :math:`\rho = \infty,\; \alpha = 1,\; z_{\rm seed}=0`

   do iteration

      :math:`\circ` Seed equation

      :math:`\rho^{\rm old} = \rho,\; \rho = {\tilde {\boldsymbol r}}^* \cdot {\boldsymbol r}`

      :math:`\beta = \rho / \rho^{\rm old}`

      :math:`{\boldsymbol q} = (z_{\rm seed} {\hat I} - {\hat H}){\boldsymbol r}`

      :math:`\alpha^{\rm old} = \alpha,\; \alpha = \frac{\rho}{{\tilde {\boldsymbol r}}^*\cdot{\boldsymbol q} - \beta \rho / \alpha }`

      :math:`\circ` Shifted equation

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

      Search :math:`k` which gives the smallest :math:`|\pi_k|` . :math:`\rightarrow z_{\rm seed},\; \pi_{\rm seed},\; \pi_{\rm seed}^{\rm old}`

      :math:`{\boldsymbol r} = {\boldsymbol r} / \pi_{\rm seed},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r}^{\rm old} / \pi_{\rm seed}^{\rm old},\; {\tilde {\boldsymbol r}} = {\tilde {\boldsymbol r}} / \pi_{\rm seed}^*,\; {\tilde {\boldsymbol r}}^{\rm old} = {\tilde {\boldsymbol r}}^{\rm old} / \pi_{\rm seed}^{{\rm old}*}`

      :math:`\alpha = (\pi_{\rm seed}^{\rm old} / \pi_{\rm seed}) \alpha`, :math:`\rho = \rho / (\pi_{\rm seed}^{\rm old} \pi_{\rm seed}^{\rm old})`

      :math:`\{\pi_k = \pi_k / \pi_{\rm seed},\; \pi_k^{\rm old} = \pi_k^{\rm old} / \pi_{\rm seed}^{\rm old}\}`

      if( :math:`|{\boldsymbol r}| <` Threshold) exit

   end do iteration

end do :math:`j`

Shifted COCG method with seed switching technique
-------------------------------------------------

This method is obtained by
:math:`{\tilde {\boldsymbol r}} = {\boldsymbol r}^*,\; {\tilde {\boldsymbol r}}^{\rm old} = {\boldsymbol r}^{{\rm old}*}`
in the BiCG method.

:math:`G_{i j}(z_k) = 0 (i=1 \cdots N_L,\; j = 1 \cdots N_R,\; k=1 \cdots N_z)`

do :math:`j = 1 \cdots N_R`

   :math:`{\boldsymbol r} = {\boldsymbol \varphi_j}`, :math:`{\boldsymbol r}^{\rm old} = {\bf 0}`

   :math:`p_{i k} = 0(i=1 \cdots N_L,\; k=1 \cdots N_z),\; \pi_k=\pi_k^{\rm old} = 1(k=1 \cdots N_z)`

   :math:`\rho = \infty,\; \alpha = 1,\; z_{\rm seed}=0`

   do iteration

      :math:`\circ` Seed equation

      :math:`\rho^{\rm old} = \rho,\; \rho = {\boldsymbol r} \cdot {\boldsymbol r}`

      :math:`\beta = \rho / \rho^{\rm old}`

      :math:`{\boldsymbol q} = (z_{\rm seed} {\hat I} - {\hat H}){\boldsymbol r}`

      :math:`\alpha^{\rm old} = \alpha,\; \alpha = \frac{\rho}{{\boldsymbol r}\cdot{\boldsymbol q} - \beta \rho / \alpha }`

      :math:`\circ` Shifted equation

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

      Search :math:`k` which gives the smallest :math:`|\pi_k|` . :math:`\rightarrow z_{\rm seed},\; \pi_{\rm seed},\; \pi_{\rm seed}^{\rm old}`
                  
      :math:`{\boldsymbol r} = {\boldsymbol r} / \pi_{\rm seed},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r}^{\rm old} / \pi_{\rm seed}^{\rm old}`

      :math:`\alpha = (\pi_{\rm seed}^{\rm old} / \pi_{\rm seed}) \alpha`, :math:`\rho = \rho / (\pi_{\rm seed}^{\rm old} \pi_{\rm seed}^{\rm old})`

      :math:`\{\pi_k = \pi_k/\pi_{\rm seed},\; \pi_k^{\rm old} = \pi_k^{\rm old} / \pi_{\rm seed}^{\rm old}\}`

      if( :math:`|{\boldsymbol r}| <` Threshold) exit

   end do iteration

end do :math:`j`

Shifted CG method with seed switching technique
-----------------------------------------------

This method is obtained by
:math:`{\tilde {\boldsymbol r}} = {\boldsymbol r},\; {\tilde {\boldsymbol r}}^{\rm old} = {\boldsymbol r}^{\rm old}`
in the BiCG method.

:math:`G_{i j}(z_k) = 0 (i=1 \cdots N_L,\; j = 1 \cdots N_R,\; k=1 \cdots N_z)`

do :math:`j = 1 \cdots N_R`

   :math:`{\boldsymbol r} = {\boldsymbol \varphi_j}`, :math:`{\boldsymbol r}^{\rm old} = {\bf 0}`

   :math:`p_{i k} = 0(i=1 \cdots N_L,\; k=1 \cdots N_z),\; \pi_k=\pi_k^{\rm old} = 1(k=1 \cdots N_z)`

   :math:`\rho = \infty,\; \alpha = 1,\; z_{\rm seed}=0`

   do iteration

      :math:`\circ` Seed equation

      :math:`\rho^{\rm old} = \rho,\; \rho = {\boldsymbol r}^* \cdot {\boldsymbol r}`

      :math:`\beta = \rho / \rho^{\rm old}`

      :math:`{\boldsymbol q} = (z_{\rm seed} {\hat I} - {\hat H}){\boldsymbol r}`

      :math:`\alpha^{\rm old} = \alpha,\; \alpha = \frac{\rho}{{\boldsymbol r}^* \cdot {\boldsymbol q} - \beta \rho / \alpha }`

      :math:`\circ` Shifted equation

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

      Search :math:`k` which gives the minimum value of :math:`|\pi_k|` . :math:`\rightarrow z_{\rm seed},\; \pi_{\rm seed},\; \pi_{\rm seed}^{\rm old}`

      :math:`{\boldsymbol r} = {\boldsymbol r} / \pi_{\rm seed},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r}^{\rm old} / \pi_{\rm seed}^{\rm old}`

      :math:`\alpha = (\pi_{\rm seed}^{\rm old} / \pi_{\rm seed}) \alpha`, :math:`\rho = \rho / {\pi_{\rm seed}^{\rm old}}^2`

      :math:`\{\pi_k = \pi_k/\pi_{\rm seed},\; \pi_k^{\rm old} = \pi_k^{\rm old}/\pi_{\rm seed}^{\rm old}\}`

      if( :math:`|{\boldsymbol r}| <` Threshold) exit

   end do iteration

end do :math:`j`

