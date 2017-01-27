Usage
=====

The calculation is done to utilize functions by the following
procedures.

-  Initialization (:ref:`init`)

-  Update results iteratively (:ref:`update`)

-  (Optional) Take the information for the restart (:ref:`getcoef`, :ref:`getvec`)

-  Finalization (:ref:`finalize`)

The restart calculation can be done by the following procedures.

-  Initialization with the information of the previous calculation (:ref:`restart`)

-  Update results iteratively (:ref:`update`)

-  (Optional) Take the information for the further restart (:ref:`getcoef`, :ref:`getvec`)

-  Finalization (:ref:`finalize`)

.. warning::

   Since :math:`K\omega` is **not** thread safe,
   these routine must be called from the outside of the OpenMP-parallel region.
   
For FORTRAN, the modules can be called by

.. code-block:: fortran

   USE komega_cg_r ! Conjugate-gradient method for real vectors
   USE komega_cg_c ! Conjugate-gradient method for complex vectors
   USE komega_cocg ! Conjugate-orthogonal conjugate-gradient mehod
   USE komega_bicg ! Biconjugate-gradient method

To utilize routines of
MPI / Hybrid parallelization version, the modules can be called as folows:

.. code-block:: fortran

   USE pkomega_cg_r
   USE pkomega_cg_c
   USE pkomega_cocg
   USE pkomega_bicg

When we call :math:`K\omega` from C/C++ codes,
we should include the header file as

.. code-block:: c

    #include komega_cg_r.h
    #include komega_cg_c.h
    #include komega_cocg.h
    #include komega_bicg.h

Scaler arguments should be passed as pointers.
For MPI/Hybrid parallelized routine,
the above line becomes

.. code-block:: c

    #include pkomega_cg_r.h
    #include pkomega_cg_c.h
    #include pkomega_cocg.h
    #include pkomega_bicg.h

Also the communicator argument for the routine should be
transformed from the C/C++'s one to the fortran's one as follows.

.. code-block:: c

      comm_f = MPI_Comm_c2f(comm_c);

Details of each routines
------------------------

.. _init:

\*_init
~~~~~~~

Set and initialize internal variables in libraries. These routines
should be called first before solving the shifted equation.

Syntax

   Fortran (Serial/OpenMP)

   .. code-block:: fortran

       CALL komega_cg_r_init(ndim, nl, nz, x, z, itermax, threshold)
       CALL komega_cg_c_init(ndim, nl, nz, x, z, itermax, threshold)
       CALL komega_cocg_init(ndim, nl, nz, x, z, itermax, threshold)
       CALL komega_bicg_init(ndim, nl, nz, x, z, itermax, threshold)

   Fortran (MPI/Hybrid parallel)

   .. code-block:: fortran

       CALL pkomega_cg_r_init(ndim, nl, nz, x, z, itermax, threshold, comm)
       CALL pkomega_cg_c_init(ndim, nl, nz, x, z, itermax, threshold, comm)
       CALL pkomega_cocg_init(ndim, nl, nz, x, z, itermax, threshold, comm)
       CALL pkomega_bicg_init(ndim, nl, nz, x, z, itermax, threshold, comm)

   C/C++ Serial/OpenMP

   .. code-block:: c

       komega_cg_r_init(&ndim, &nl, &nz, x, z, &itermax, &threshold);
       komega_cg_c_init(&ndim, &nl, &nz, x, z, &itermax, &threshold);
       komega_cocg_init(&ndim, &nl, &nz, x, z, &itermax, &threshold);
       komega_bicg_init(&ndim, &nl, &nz, x, z, &itermax, &threshold);

   C/C++ MPI/Hybrid parallel

   .. code-block:: c

       pkomega_cg_r_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);
       pkomega_cg_c_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);
       pkomega_cocg_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);
       pkomega_bicg_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);

Parameters

   .. code-block:: fortran

      INTEGER,INTENT(IN) :: ndim
   ..

      The dimension of solution vectors for the linearized equation.
      ``ndim`` for the dimension of variables in other routine is
      equal to this.

   .. code-block:: fortran

      INTEGER,INTENT(IN) :: nl
   ..

      The dimension of projected solution vectors.
      ``nl`` for the dimension of variables in other routine is
      equal to this.

   .. code-block:: fortran
                
      INTEGER,INTENT(IN) :: nz
   ..

      The number of shifted points.
      ``nz`` for the dimension of variables in other routine is
      equal to this.

   .. code-block:: fortran

      REAL(8),INTENT(OUT) :: x(nl*nz) ! (for "CG_R_init", "CG_C_init")
      COMPLEX(8),INTENT(OUT) :: x(nl*nz) ! (for other cases)
   ..

      The solution vector. In this procedure, ``0`` vector is returned.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: z(nz) ! (for "CG_R_init", "CG_C_init")
      COMPLEX(8),INTENT(IN) :: z(nz) ! (for other cases)
   ..

      Shifted points.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: itermax
   ..

      The maximum iteration number for allocating arrays for the restart calculation.
      When ``itermax=0`` , these arrays are not allocated,
      and the restart calculation described later becomes unavailable.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: threshold
   ..

      The threshold value for the convergence determination.
      When the 2-norm of the residual vector for the seed equation
      becomes smaller than this value, the calculation is finished.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: comm
   ..

      Only for MPI / Hybrid parallelization
      version. Communicators for MPI such as ``MPI_COMM_WORLD``.

.. _restart:
   
\*_restart
~~~~~~~~~~

For the restart calculation, these routines are used instead of :ref:`init`.
Set and initialize internal variables in libraries.
These routines should be called first before solving the shifted equation.

Syntax

   Fortran (Serial/OpenMP)

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

   Fortran (MPI/hybrid parallel)

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

   C/C++ (Serial/OpenMP)

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

   C/C++ (MPI/hybrid parallel)

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

Parameters

   .. code-block:: fortran

      INTEGER,INTENT(IN) :: ndim
      INTEGER,INTENT(IN) :: nl
      INTEGER,INTENT(IN) :: nz
      REAL(8),INTENT(OUT) :: x(nl*nz)
      REAL(8),INTENT(IN) :: z(nz) ! (for "CG_R_restart", "CG_C_restart")
      COMPLEX(8),INTENT(IN) :: z(nz) ! (Other)
      INTEGER,INTENT(IN) :: itermax
      REAL(8),INTENT(IN) :: threshold
      INTEGER,INTENT(IN) :: comm
   ..
   
      The definition is same as :ref:`init`. See the parameters in :ref:`init`.

   .. code-block:: fortran

      INTEGER,INTENT(OUT) :: status(3)
   ..
   
      The error code is returned.

      First component(``status(1)``)
      
         If the solution is converged or a breakdown occurs,
         the current total number of iteration with the minus sign is returned.
         In other cases, this routine returns the current total number of iteration.
         The calculation is continuable only when ``status(1)`` is the positive value;
         otherwise the result is meaningless even if the calculation is continued.

      Second component(``status(2)``)
      
         ``1`` is returned if ``itermax`` is set as a finite value and the
         convergence condition is not satisfied at the ``itermax``\ -th iteration.
         ``2`` is returned if :math:`\alpha` diverges.
         ``3`` is returned if :math:`\pi_{\rm seed}` becomes 0.
         In the case of ``COCG_restart`` or ``BiCG_restart``,
         ``4`` is returned if the residual vector and the shadow residual vector are orthogonal.
         In other cases, ``0`` is returned.

      Third component(``status(3)``)
      
         The index of the seed point is returned.

   .. code-block:: fortran
                   
      INTEGER,INTENT(IN) :: iter_old
   ..
   
      The number of iteration for the previous calculation.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: v2(ndim) ! (for "CG_R_restart")
      COMPLEX(8),INTENT(IN) :: v2(ndim) ! (Other)
   ..
   
      The residual vector at the last step for the previous calculation.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: v12(ndim) ! (for "CG_R_restart")
      COMPLEX(8),INTENT(IN) :: v12(ndim) ! (Other)
   ..

      The residual vector at the second from the last step for the previous calculation.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: alpha_save(iter_old) ! (for "CG_R_restart", "CG_C_restart")
      COMPLEX(8),INTENT(IN) :: alpha_save(iter_old) ! (Other)
   ..                   

      The parameters :math:`\alpha` obtained by the
      previous calculation at each steps by (Bi)CG methods.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: beta_save(iter_old) ! (for "CG_R_restart", "CG_C_restart")
      COMPLEX(8),INTENT(IN) :: beta_save(iter_old) ! (Other)
   ..                   

      The parameters :math:`\beta` obtained
      by the previous calculation at each steps by (Bi)CG methods.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: z_seed ! (for "CG_R_restart", "CG_C_restart")
      COMPLEX(8),INTENT(IN) :: z_seed ! (Other)
   ..                   

      The value of the seed shift for the previous calculation.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: r_l_save(nl,iter_old) ! (for "CG_R_restart")
      COMPLEX(8),INTENT(IN) :: r_l_save(nl,iter_old) ! (Other)
   ..                   

      The projected residual vector at each iteration for the previous calculation.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: v4(ndim) ! (for "CG_R_restart")
      COMPLEX(8),INTENT(IN) :: v4(ndim) ! (Other)
   ..
   
      The shadow residual vector at the last step for the previous calculation.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: v14(ndim) ! (for "CG_R_restart")
      COMPLEX(8),INTENT(IN) :: v14(ndim) ! (Other)
   ..

      The shadow residual vector at the second last step for the previous calculation.

.. _update:
      
\*_update
~~~~~~~~~

It is called alternately with the matrix-vector product
in the loop and updates the solution.

Syntax

   Fortran (Serial/OpenMPI)

   .. code-block:: fortran

       CALL komega_cg_r_update(v12, v2, x, r_l, status)
       CALL komega_cg_c_update(v12, v2, x, r_l, status)
       CALL komega_cocg_update(v12, v2, x, r_l, status)
       CALL komega_bicg_update(v12, v2, v14, v4, x, r_l, status)

   Fortran (MPI/hybrid parallel)

   .. code-block:: fortran

       CALL pkomega_cg_r_update(v12, v2, x, r_l, status)
       CALL pkomega_cg_c_update(v12, v2, x, r_l, status)
       CALL pkomega_cocg_update(v12, v2, x, r_l, status)
       CALL pkomega_bicg_update(v12, v2, v14, v4, x, r_l, status)

   C/C++ (Serial/OpenMPI)

   .. code-block:: c

       komega_cg_r_update(v12, v2, x, r_l, status);
       komega_cg_c_update(v12, v2, x, r_l, status);
       komega_cocg_update(v12, v2, x, r_l, status);
       komega_bicg_update(v12, v2, v14, v4, x, r_l, status);

   C/C++ (MPI/hybrid parallel)

   .. code-block:: c

       pkomega_cg_r_update(v12, v2, x, r_l, status);
       pkomega_cg_c_update(v12, v2, x, r_l, status);
       pkomega_cocg_update(v12, v2, x, r_l, status);
       pkomega_bicg_update(v12, v2, v14, v4, x, r_l, status);

Parameters

   .. code-block:: fortran

      REAL(8),INTENT(INOUT) :: v12(ndim) ! (for "CG_R_update")
      COMPLEX(8),INTENT(INOUT) :: v12(ndim) ! (Other)
   ..

      The product of the residual vector (``v2``) and the matrix.
      This routine returns the 2-norm of the updated residual vector
      as a first element of this array.
      This returned value is used, for examples, for printing the convergence profile.

   .. code-block:: fortran

      REAL(8),INTENT(INOUT) :: v2(ndim) ! (for "CG_R_update")
      COMPLEX(8),INTENT(INOUT) :: v2(ndim) ! (Other)
   ..
   
      The residual vector is input and the updated residual vector is output.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: v14(ndim) ! (for "CG_R_update")
      COMPLEX(8),INTENT(IN) :: v14(ndim) ! (Other)
   ..

      The product of the shadow residual vector (``v4``) and the matrix is input.

   .. code-block:: fortran

      REAL(8),INTENT(INOUT) :: v4(ndim) ! (for "CG_R_update")
      COMPLEX(8),INTENT(INOUT) :: v4(ndim) ! (Other)
   ..

      The shadow residual vector is input and the updated vector is output.

   .. code-block:: fortran

      INTEGER,INTENT(OUT) :: status(3)
   ..
   
      The error code is returned.

      First component (``status(1)``)
      
         If the solution is converged or a breakdown occurs,
         the current total number of iteration with the minus sign is returned.
         In other cases,
         this routine returns the current total number of iteration.
         The calculation is continuable only when ``status(1)`` is the positive value;
         otherwise the result is meaningless even if the calculation is continued.

      Second component (``status(2)``)
      
         ``1`` is returned if ``itermax`` is set as a finite value in the
         :ref:`init` routine and the convergence condition is not satisfied
         at the ``itermax``\ -th iteration.
         ``2`` is returned if :math:`\alpha` diverges.
         ``3`` is returned if :math:`\pi_{\rm seed}` becomes 0.
         In the case of ``COCG_update`` or ``BiCG_update``,
         ``4`` is returned if the residual vector and
         the shadow residual vector are orthogonal.
         In other cases, ``0`` is returned.

      Third component (``status(3)``)
   
         The index of the seed point is returned.

.. _getcoef:
         
\*_getcoef
~~~~~~~~~~

Get the coefficients used in the restart calculation.
To call these routines,
``itermax`` in :ref:`init` routine must not be ``0`` .

The total number of iteration (``iter_old``) used in this routine
is computed by using ``status`` which is an output of :ref:`update` as follows:

.. code-block:: fortran

   iter_old = ABS(status(1))

Syntax

   Fortran (Serial/OpenMP)

   .. code-block:: fortran

       CALL komega_cg_r_getcoef(alpha_save, beta_save, z_seed, r_l_save)
       CALL komega_cg_c_getcoef(alpha_save, beta_save, z_seed, r_l_save)
       CALL komega_cocg_getcoef(alpha_save, beta_save, z_seed, r_l_save)
       CALL komega_bicg_getcoef(alpha_save, beta_save, z_seed, r_l_save)

   Fortran (MPI/hybrid parallel)

   .. code-block:: fortran

       CALL pkomega_cg_r_getcoef(alpha_save, beta_save, z_seed, r_l_save)
       CALL pkomega_cg_c_getcoef(alpha_save, beta_save, z_seed, r_l_save)
       CALL pkomega_cocg_getcoef(alpha_save, beta_save, z_seed, r_l_save)
       CALL pkomega_bicg_getcoef(alpha_save, beta_save, z_seed, r_l_save)

   C/C++ (Serial/OpenMP)

   .. code-block:: c

       komega_cg_r_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
       komega_cg_c_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
       komega_cocg_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
       komega_bicg_getcoef(alpha_save, beta_save, &z_seed, r_l_save);

   C/C++ (MPI/hybrid parallel)

   .. code-block:: c

       pkomega_cg_r_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
       pkomega_cg_c_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
       pkomega_cocg_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
       pkomega_bicg_getcoef(alpha_save, beta_save, &z_seed, r_l_save);

Parameters

   .. code-block:: fortran

      REAL(8),INTENT(OUT) :: alpha_save(iter_old) ! (for "CG_R_restart", "CG_C_restart")
      COMPLEX(8),INTENT(OUT) :: alpha_save(iter_old) ! (Other)
   ..
   
      The parameters :math:`\alpha` of the (Bi)CG method at each iteration.

   .. code-block:: fortran

      REAL(8),INTENT(OUT) :: beta_save(iter_old) ! (for "CG_R_restart", "CG_C_restart")
      COMPLEX(8),INTENT(OUT) :: beta_save(iter_old) ! (Other)
   ..                   

      The parameters :math:`\beta` of the (Bi)CG method at each iteration.

   .. code-block:: fortran

      REAL(8),INTENT(OUT) :: z_seed ! (for "CG_R_restart", "CG_C_restart")
      COMPLEX(8),INTENT(OUT) :: z_seed ! (Other)
   ..                   

      Seed shift.

   .. code-block:: fortran

      REAL(8),INTENT(IN) :: r_l_save(nl,iter_old) ! ("CG_R_restart")
      COMPLEX(8),INTENT(IN) :: r_l_save(nl,iter_old) ! (Other)
   ..                   

      The projected residual vectors at each iteration.

.. _getvec:
      
\*_getvec
~~~~~~~~~

Get the residual vectors to use the restart calculation.
To call these routines,
``itermax`` in the :ref:`init` routine must not be ``0``.

Syntax

   Fortran (Serial/OpenMP)

   .. code-block:: fortran

       CALL komega_cg_r_getvec(r_old)
       CALL komega_cg_c_getvec(r_old)
       CALL komega_cocg_getvec(r_old)
       CALL komega_bicg_getvec(r_old, r_tilde_old)

   Fortran (MPI/hybrid parallel)

   .. code-block:: fortran

       CALL pkomega_cg_r_getvec(r_old)
       CALL pkomega_cg_c_getvec(r_old)
       CALL pkomega_cocg_getvec(r_old)
       CALL pkomega_bicg_getvec(r_old, r_tilde_old)

   C/C++ (Serial/OpenMP)

   .. code-block:: c

       komega_cg_r_getvec(r_old);
       komega_cg_c_getvec(r_old);
       komega_cocg_getvec(r_old);
       komega_bicg_getvec(r_old, r_tilde_old);

   C/C++ (MPI/hybrid parallel)

   .. code-block:: c

       pkomega_cg_r_getvec(r_old);
       pkomega_cg_c_getvec(r_old);
       pkomega_cocg_getvec(r_old);
       pkomega_bicg_getvec(r_old, r_tilde_old);

Parameters

   .. code-block:: fortran

      REAL(8),INTENT(OUT) :: r_old(ndim) ! (for "CG_R_getvec")
      COMPLEX(8),INTENT(OUT) :: r_old(ndim) ! (Other)
   ..

      The residual vector at the second last step in the previous calculation.

   .. code-block:: fortran

      COMPLEX(8),INTENT(OUT) :: r_tilde_old(ndim)
   ..

      The shadow residual vector at the second last step in the previous calculation.

\*_getresidual
~~~~~~~~~~~~~~

Get the values of 2-norm of the residual vector at each shift points.
These routines can be called from anywhere between :ref:`init`
and :ref:`finalize` .
These routines do not affect the calculation results.

Syntax

   Fortran (Serial/OpenMP)

   .. code-block:: fortran

       CALL komega_cg_r_getresidual(res)
       CALL komega_cg_c_getresidual(res)
       CALL komega_cocg_getresidual(res)
       CALL komega_bicg_getresidual(res)

   Fortran (MPI/hybrid parallel)

   .. code-block:: fortran

       CALL pkomega_cg_r_getresidual(res)
       CALL pkomega_cg_c_getresidual(res)
       CALL pkomega_cocg_getresidual(res)
       CALL pkomega_bicg_getresidual(res)

   C/C++ (Serial/OpenMP)

   .. code-block:: c

       komega_cg_r_getresidual(res);
       komega_cg_c_getresidual(res);
       komega_cocg_getresidual(res);
       komega_bicg_getresidual(res);

   C/C++ (MPI/hybrid parallel)

   .. code-block:: c

       pkomega_cg_r_getresidual(res);
       pkomega_cg_c_getresidual(res);
       pkomega_cocg_getresidual(res);
       pkomega_bicg_getresidual(res);

Parameters

   .. code-block:: fortran

      COMPLEX(8),INTENT(OUT) :: res(nz)
   ..

      The values of 2-norm of the residual vector at each shift points are
      returned.

.. _finalize:
      
\*_finalize
~~~~~~~~~~~

Release memories of the arrays stored in the library.

Syntax

   Fortran (Serial/OpenMP)

   .. code-block:: fortran

       CALL komega_cg_r_finalize()
       CALL komega_cg_c_finalize()
       CALL komega_cocg_finalize()
       CALL komega_bicg_finalize()

   Fortran (MPI/hybrid parallel)

   .. code-block:: fortran

       CALL pkomega_cg_r_finalize()
       CALL pkomega_cg_c_finalize()
       CALL pkomega_cocg_finalize()
       CALL pkomega_bicg_finalize()

   C/C++ (Serial/OpenMP)

   .. code-block:: c

       komega_cg_r_finalize();
       komega_cg_c_finalize();
       komega_cocg_finalize();
       komega_bicg_finalize();

   C/C++ (MPI/hybrid parallel)

   .. code-block:: c

       pkomega_cg_r_finalize();
       pkomega_cg_c_finalize();
       pkomega_cocg_finalize();
       pkomega_bicg_finalize();

Sample codes for using shifted BiCG library
-------------------------------------------

As a typical example, the usage of shifted BiCG library is shown below.

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

