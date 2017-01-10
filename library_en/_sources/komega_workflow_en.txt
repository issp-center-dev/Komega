Schematic workflow of this library
==================================

In the following description, the loop for :math:`N_R` is omitted for
simplicity and instead of :math:`G_{i j}(z_k)`,
the :math:`N_L`\ -dimensional vector :math:`{\bf x}_{k}`
is obtained by using the library.

The essential procedures to use the library.

-  Give the vector size :math:`N_H` corresponding to the size of the
   Hilbert space and the number of the frequency :math:`z`.

-  Allocate the two vectors (in the case of BiCG method, four vectors)
   with the size of :math:`N_H`.

-  Give the function for the Hamiltonian-vector production.

-  Allocate the solution vectors. It is noted that the number of the
   solution vectors is not always equal to :math:`N_H \times N_z`. In
   fact, the number in the previous section is :math:`N_L \times N_z`.
   In this case, the (bi-)conjugate gradient vector :math:`{{\bf p}_k}`
   is the vector with :math:`N_L`-dimension and :math:`N_z` array. The
   transformation of the dimension of the residual vector from
   :math:`N_H` to :math:`N_L` must be done explicitly.

   .. math::

      \begin{aligned}
          {\bf r}^{\rm L} = {\hat P}^\dagger {\boldsymbol r}, \qquad
          {\hat P} \equiv ({\boldsymbol \varphi}_1, \cdots, {\boldsymbol \varphi}_{N_L})
        \end{aligned}

-  If the result converges, ``komega_????_update`` return the first element of ``status``
   as a negative integer.
   Therefore, please exit loop when ``status(1) < 0`` .

-  The 2-norm is used for the convergence check in the routine ``komega_????_update``.
   Therefore, if 2-norms of residual vectors at all shift points becomes smaller than K,
   this routine assumes the result is converged.

-  Conserve :math:`\alpha, \beta, {\bf r}^{\rm L}` at each iteration and
   to use the above vectors later, set the maximum number of iteration
   ``itermax``.

The names of the routines is defined as follows.

-  ``BiCG_init``, ``COCG_init``, ``CG_C_init``, ``CG_R_init``

   Set the initial conditions such as the allocation of variables used
   in the library.

-  ``BiCG_update``, ``COCG_update``, ``CG_C_update``, ``CG_R_update``

   These routines are called in the iteration to update the solution
   vectors.

-  ``BiCG_finalize``, ``COCG_finalize``, ``CG_C_finalize``,
   ``CG_R_finalize``

   Release the allocated vectors in the library.

-  ``BiCG_getcoef``, ``COCG_getcoef``, ``CG_C_getcoef``,
   ``CG_R_getcoef``

   Get the :math:`\alpha`, :math:`\beta`, :math:`z_{\rm seed}`,
   :math:`{\bf r}^{\rm L}` conserved at each iteration.

-  ``BiCG_getvec``, ``COCG_getvec``, ``CG_C_getvec``, ``CG_R_getvec``

   Get the vectors :math:`{\boldsymbol r}`,
   :math:`{\boldsymbol r}^{\rm old}`, :math:`{\tilde {\boldsymbol r}}`,
   :math:`{\tilde {\boldsymbol r}}^{\rm old}`.

-  ``BiCG_restart``, ``COCG_restart``, ``CG_C_restart``,
   ``CG_R_restart``

The schematic workflow of shifted BiCG library
----------------------------------------------

Allocate :math:`{\boldsymbol v}_{1 2}`, :math:`{\boldsymbol v}_{1 3}`,
:math:`{\boldsymbol v}_2`, :math:`{\boldsymbol v}_3`,
:math:`\{{\bf x}_k\}, {\bf r}^{\rm L}`
:math:`{\boldsymbol v}_2 = {\boldsymbol \varphi_j}`

``komega_BiCG_init(N_H, N_L, N_z, x, z, itermax, threshold)`` start

   Allocate :math:`{\boldsymbol v}_3`, :math:`{\boldsymbol v}_5`,
   :math:`\{\pi_k\}` , :math:`\{\pi_k^{\rm old}\}`, :math:`\{{\bf p}_k\}`

   Copy :math:`\{z_k\}`

   If ``itermax`` :math:`\neq` ``0`` ,
   allocate arrays to store :math:`\alpha`, :math:`\beta`,
   and:math:`{\bf r}^{\rm L}` at each iteration.

   :math:`{\boldsymbol v}_4 = {\boldsymbol v}_2^*` (an arbitrary vector),
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
   [Or :math:`({\boldsymbol v}_{1 2}, {\boldsymbol v}_{1 4}) = {\hat H} ({\boldsymbol v}_2, {\boldsymbol v}_4)` ]

   ``komega_BiCG_update(v_12, v_2, v_14, v_4, x, r_small, status)`` start

      :math:`\circ` Seed equation

      :math:`\rho^{\rm old} = \rho,\; \rho = {\boldsymbol v}_4^* \cdot {\boldsymbol v}_2`

      :math:`\beta = \rho / \rho^{\rm old}`

      :math:`{\boldsymbol v}_{1 2} = z_{\rm seed} {\boldsymbol v}_2 - {\boldsymbol v}_{1 2}`,
      :math:`{\boldsymbol v}_{1 4} = z_{\rm seed}^* {\boldsymbol v}_4 - {\boldsymbol v}_{1 4}`

      :math:`\alpha^{\rm old} = \alpha,\; \alpha = \frac{\rho}{{\boldsymbol v}_3^* \cdot {\boldsymbol v}_{1 2} - \beta \rho / \alpha }`

      :math:`\circ` Shifted equation

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

      Search :math:`k` which gives the smallest :math:`|\pi_k|` .
      :math:`\rightarrow z_{\rm seed},\; \pi_{\rm seed},\; \pi_{\rm seed}^{\rm old}`

      :math:`{\boldsymbol v}_2 = {\boldsymbol v}_2 / \pi_{\rm seed}`,
      :math:`{\boldsymbol v}_3 = {\boldsymbol v}_3 / \pi_{\rm seed}^{\rm old}`,
      :math:`{\boldsymbol v}_4 = {\boldsymbol v}_4 / \pi_{\rm seed}^{*}`,
      :math:`{\boldsymbol v}_5 = {\boldsymbol v}_5 / \pi_{\rm seed}^{\rm old *}`

      :math:`\alpha = (\pi_{\rm seed}^{\rm old} / \pi_{\rm seed}) \alpha`,
      :math:`\rho = \rho / (\pi_{\rm seed}^{\rm old} \pi_{\rm seed}^{\rm old})`

      :math:`\{\pi_k = \pi_k / \pi_{\rm seed},\; \pi_k^{\rm old} = \pi_k^{\rm old} / \pi_{\rm seed}^{\rm old}\}`

   ``komega_BiCG_update`` finish

   if(status(1) < 0 (This indicates :math:`|{\boldsymbol v}_2| <` Threshold)) exit

end do iteration

``komega_BiCG_finalize`` start

   Deallocate :math:`{\boldsymbol v}_4`, :math:`{\boldsymbol v}_5`,
   :math:`\{\pi_k\}`, :math:`\{\pi_k^{\rm old}\}`, :math:`\{{\bf p}_k\}`

``komega_BiCG_finalize`` finish

The schematic workflow of shifted COCG library
----------------------------------------------

Allocate :math:`{\boldsymbol v}_1`, :math:`{\boldsymbol v}_2`,
:math:`\{{\bf x}_k\}, {\bf r}^{\rm L}`
:math:`{\boldsymbol v}_2 = {\boldsymbol \varphi_j}`

``COCG_init(N_H, N_L, N_z, x, z, itermax, threshold)`` start

   Allocate :math:`{\boldsymbol v}_3`, :math:`\{\pi_k\}`,
   :math:`\{\pi_k^{\rm old}\}`, :math:`\{{\bf p}_k\}`

   Copy :math:`\{z_k\}`

   If ``itermax`` :math:`\neq` ``0`` , allocate arrays
   to store :math:`\alpha`, :math:`\beta`, and :math:`{\bf r}^{\rm L}` .

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

      :math:`\circ` Seed equationw

      :math:`\rho^{\rm old} = \rho,\; \rho = {\boldsymbol v}_2 \cdot {\boldsymbol v}_2`

      :math:`\beta = \rho / \rho^{\rm old}`

      :math:`{\boldsymbol v}_1 = z_{\rm seed} {\boldsymbol v}_2 - {\boldsymbol v}_1`

      :math:`\alpha^{\rm old} = \alpha,\; \alpha = \frac{\rho}{{\boldsymbol v}_2 \cdot {\boldsymbol v}_1 - \beta \rho / \alpha }`

      :math:`\circ` Shifted equations

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

      Search :math:`k` which gives the smallest `|\pi_k|` .
      :math:`\rightarrow z_{\rm seed},\; \pi_{\rm seed},\; \pi_{\rm seed}^{\rm old}`

      :math:`{\boldsymbol v}_2 = {\boldsymbol v}_2 / \pi_{\rm seed}`,
      :math:`{\boldsymbol v}_3 = {\boldsymbol v}_3 / \pi_{\rm seed}^{\rm old}`

      :math:`\alpha = (\pi_{\rm seed}^{\rm old} / \pi_{\rm seed}) \alpha`,
      :math:`\rho = \rho / (\pi_{\rm seed}^{\rm old} \pi_{\rm seed}^{\rm old})`

      :math:`\{\pi_k = \pi_k / \pi_{\rm seed},\; \pi_k^{\rm old} = \pi_k^{\rm old} / \pi_{\rm seed}^{\rm old}\}`

   ``COCG_update`` finish

   if(status(1) < 0 (This indicates :math:`|{\boldsymbol v}_2| <` Threshold.)) exit

end do iteration

``COCG_finalize`` start

   Deallocate :math:`{\boldsymbol v}_3`, :math:`\{\pi_k\}`,
   :math:`\{\pi_k^{\rm old}\}`, :math:`\{{\bf p}_k\}`

``COCG_finalize`` finish

The schematic workflow of shifted CG library
--------------------------------------------

The workflow is the same as that of the shifted COCG library.
