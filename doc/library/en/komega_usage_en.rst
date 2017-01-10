Usage
=====

The calculation is done to utilize functions by the following
procedures.

-  Initialization (init function)

-  Update (update function)

-  Output numerical results (call getcoef, getvec functions and output
   informations)

-  Finalization (finalize function)

The restart calculation can be done by the following procedures.

-  Initialization(restart function)

-  Update (update function)

-  Output numerical results (call getcoef, getvec functions and output
   informations)

-  Finalization (finalize function)

For FORTRAN, the modules can be called by

.. code-block:: fortran

      USE komega_????

``"????"`` is selected from the following words(methods)
``"CG_R"``, ``"CG_C"``, ``"COCG"``, ``"BiCG"``. To utilize routines of
MPI / Hybrid parallelization version, the modules can be called as folows:

.. code-block:: fortran

      USE pkomega_????

When we call :math:`K\omega` from C/C++ codes,
we should include the header file as

.. code-block:: c

    #include komega_????.h

Scaler arguments should be passed as pointers.
For MPI/Hybrid parallelized routine,
the above line becomes

.. code-block:: c

    #include pkomega_????.h

Also the communicator argument for the routine should be
transformed from the C/C++'s one to the fortran's one as follows.

.. code-block:: c

      comm_f = MPI_Comm_c2f(comm_c);

Details of each routines
------------------------

``????_init``, ``p????_init``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set and initialize internal variables in libraries. These routines
should be called first before solving the shifted equation.

Syntax

Fortran (Serial/OpenMP)

.. code-block:: fortran

    CALL komega_CG_R_init(ndim, nl, nz, x, z, itermax, threshold)
    CALL komega_CG_C_init(ndim, nl, nz, x, z, itermax, threshold)
    CALL komega_COCG_init(ndim, nl, nz, x, z, itermax, threshold)
    CALL komega_BiCG_init(ndim, nl, nz, x, z, itermax, threshold)

Fortran (MPI/Hybrid parallel)

.. code-block:: fortran

    CALL pkomega_CG_R_init(ndim, nl, nz, x, z, itermax, threshold, comm)
    CALL pkomega_CG_C_init(ndim, nl, nz, x, z, itermax, threshold, comm)
    CALL pkomega_COCG_init(ndim, nl, nz, x, z, itermax, threshold, comm)
    CALL pkomega_BiCG_init(ndim, nl, nz, x, z, itermax, threshold, comm)

C/C++ Serial/OpenMP

.. code-block:: c

    komega_CG_R_init(&ndim, &nl, &nz, x, z, &itermax, &threshold);
    komega_CG_C_init(&ndim, &nl, &nz, x, z, &itermax, &threshold);
    komega_COCG_init(&ndim, &nl, &nz, x, z, &itermax, &threshold);
    komega_BiCG_init(&ndim, &nl, &nz, x, z, &itermax, &threshold);

C/C++ MPI/Hybrid parallel

.. code-block:: c

    pkomega_CG_R_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);
    pkomega_CG_C_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);
    pkomega_COCG_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);
    pkomega_BiCG_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);

Parameters

-  ``ndim``

   ``INTEGER``. Scalar. Input.
   The dimension of solution vectors for the linearized equation.

-  ``nl``

   ``INTEGER``. Scalar. Input.
   The dimension of projected solution vectors.

-  ``nz``

   ``INTEGER``. Scalar. Input. The number of shifted points.

-  ``x``

   ``DOUBLE PRECISION`` (for ``CG_R_init``), ``DOUBLE COMPLEX`` (for
   other cases). The array with the length of ``nl*nz``. Output. The
   solution vector. In this procedure, ``0`` vector is returned.

-  ``z``

   ``DOUBLE PRECISION`` (for ``CG_R_init``, ``CG_C_init``),
   ``DOUBLE COMPLEX`` (for other cases). The array with the length of
   ``nz``. Input. Shifted points.

-  ``itermax``

   ``INTEGER``. Scalar. Input.
   The maximum iteration number for allocating arrays for the restart calculation.
   When ``itermax=0`` , these arrays are not allocated,
   and the restart calculation described later becomes unavailable.

-  ``threshold``

   ``DOUBLE PRECISION``. Scalar. Input.
   The threshold value for the convergence determination.
   When the 2-norm of the residual vector for the seed equation
   becomes smaller than this value, the calculation is finished.

-  ``comm``

   ``INTEGER``. Scalar. Input. Only for MPI / Hybrid parallelization
   version. Communicators for MPI such as ``MPI_COMM_WORLD``.

``komega_????_restart``, ``pkomega_????_restart``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the restart calculation, these routines are used instead of ``?_init``.
Set and initialize internal variables in libraries.
These routines should be called first before solving the shifted equation.

Syntax

Fortran (Serial/OpenMP)

.. code-block:: fortran

    CALL komega_CG_R_restart(ndim, nl, nz, x, z, itermax, threshold, status, &
    &                 iter_old, v2, v12, alpha_save, beta_save, z_seed, r_l_save)
    CALL komega_CG_C_restart(ndim, nl, nz, x, z, itermax, threshold, status, &
    &                 iter_old, v2, v12, alpha_save, beta_save, z_seed, r_l_save)
    CALL komega_COCG_restart(ndim, nl, nz, x, z, itermax, threshold, status, &
    &                 iter_old, v2, v12, alpha_save, beta_save, z_seed, r_l_save)
    CALL komega_BiCG_restart(ndim, nl, nz, x, z, itermax, threshold, status, &
    &                 iter_old, v2, v12, v4, v14, alpha_save, beta_save, &
    &                 z_seed, r_l_save)

Fortran (MPI/hybrid parallel)

.. code-block:: fortran

    CALL pkomega_CG_R_restart(ndim, nl, nz, x, z, itermax, threshold, comm, status, &
    &                 iter_old, v2, v12, alpha_save, beta_save, z_seed, r_l_save)
    CALL pkomega_CG_C_restart(ndim, nl, nz, x, z, itermax, threshold, comm, status, &
    &                 iter_old, v2, v12, alpha_save, beta_save, z_seed, r_l_save)
    CALL pkomega_COCG_restart(ndim, nl, nz, x, z, itermax, threshold, comm, status, &
    &                 iter_old, v2, v12, alpha_save, beta_save, z_seed, r_l_save)
    CALL pkomega_BiCG_restart(ndim, nl, nz, x, z, itermax, threshold, comm, status, &
    &                 iter_old, v2, v12, v4, v14, alpha_save, beta_save, &
    &                 z_seed, r_l_save)

C/C++ (Serial/OpenMP)

.. code-block:: c

    komega_CG_R_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, status, &
    &                 &iter_old, v2, v12, alpha_save, beta_save, &z_seed, r_l_save);
     komega_CG_C_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, status, &
    &                 &iter_old, v2, v12, alpha_save, beta_save, &z_seed, r_l_save);
    komega_COCG_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, status, &
    &                 &iter_old, v2, v12, alpha_save, beta_save, &z_seed, r_l_save);
    komega_BiCG_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, status, &
    &                 &iter_old, v2, v12, v4, v14, alpha_save, beta_save, &
    &                 &z_seed, r_l_save);

C/C++ (MPI/hybrid parallel)

.. code-block:: c

    pkomega_CG_R_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm, status, &
    &                 &iter_old, v2, v12, alpha_save, beta_save, &z_seed, r_l_save);
    pkomega_CG_C_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm, status, &
    &                 &iter_old, v2, v12, alpha_save, beta_save, &z_seed, r_l_save);
    pkomega_COCG_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm, status, &
    &                 &iter_old, v2, v12, alpha_save, beta_save, &z_seed, r_l_save);
    pkomega_BiCG_restart(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm, status, &
    &                 &iter_old, v2, v12, v4, v14, alpha_save, beta_save, &
    &                 &z_seed, r_l_save);

Parameters

-  ``ndim, nl, nz, x, z, itermax, threshold, comm``

   The definition is same as ``?_init``. See the parameters in ``?_init``.

-  ``status``

   ``INTEGER``. The array with the length of ``3``. Output.
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

-  ``iter_old``

   ``INTEGER``. Scalar. Input.
   The number of iteration for the previous calculation.

-  ``v2``

   ``DOUBLE PRECISION`` (for ``CG_R_restart``),
   ``DOUBLE COMPLEX`` (for other cases).
   The array with the length of ``ndim``. Input.
   The residual vector at the last step for the previous calculation.

-  ``v12``

   ``DOUBLE PRECISION`` (for ``CG_R_restart``), ``DOUBLE COMPLEX`` (for other cases).
   The array with the length of ``ndim``. Input.
   The residual vector at the second from the last step for the previous calculation.

-  ``alpha_save``

   ``DOUBLE PRECISION`` (for ``CG_R_restart``, ``CG_C_restart``),
   ``DOUBLE COMPLEX`` (for other cases).
   The array with the length of ``iter_old``. Input.
   The parameters :math:`\alpha` obtained by the
   previous calculation at each steps by (Bi)CG methods.

-  ``beta_save``

   ``DOUBLE PRECISION`` (for ``CG_R_restart``, ``CG_C_restart``),
   ``DOUBLE COMPLEX`` (for other cases).
   The array with the length of\ ``iter_old``. Input.
   The parameters :math:`\beta` obtained
   by the previous calculation at each steps by (Bi)CG methods.

-  ``z_seed``

   ``DOUBLE PRECISION`` (for ``CG_R_restart``, ``CG_C_restart``),
   ``DOUBLE COMPLEX`` (for other cases). Scalar. Input.
   The value of the seed shift for the previous calculation.

-  ``r_l_save``

   ``DOUBLE PRECISION`` (for ``CG_R_restart``), ``DOUBLE COMPLEX`` (for other cases).
   The array with the length of ``nl*iter_old``. Input.
   The projected residual vector at each iteration for the previous calculation.

-  ``v4``

   Only used for ``BiCG_restart``. ``DOUBLE COMPLEX``.
   The array with the length of ``ndim``. Input.
   The shadow residual vector at the last step for the previous calculation.

-  ``v14``

   Only used for ``BiCG_restart``. ``DOUBLE COMPLEX``.
   The array with the length of ``ndim``. Input.
   The shadow residual vector at the second last step for the previous calculation.

``komega_????_update``, ``pkomega_????_update``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is called alternately with the matrix-vector product
in the loop and updates the solution.

Syntax

Fortran (Serial/OpenMPI)

.. code-block:: fortran

    CALL komega_CG_R_update(v12, v2, x, r_l, status)
    CALL komega_CG_C_update(v12, v2, x, r_l, status)
    CALL komega_COCG_update(v12, v2, x, r_l, status)
    CALL komega_BiCG_update(v12, v2, v14, v4, x, r_l, status)

Fortran (MPI/hybrid parallel)

.. code-block:: fortran

    CALL pkomega_CG_R_update(v12, v2, x, r_l, status)
    CALL pkomega_CG_C_update(v12, v2, x, r_l, status)
    CALL pkomega_COCG_update(v12, v2, x, r_l, status)
    CALL pkomega_BiCG_update(v12, v2, v14, v4, x, r_l, status)

C/C++ (Serial/OpenMPI)

.. code-block:: c

    komega_CG_R_update(v12, v2, x, r_l, status);
    komega_CG_C_update(v12, v2, x, r_l, status);
    komega_COCG_update(v12, v2, x, r_l, status);
    komega_BiCG_update(v12, v2, v14, v4, x, r_l, status);

C/C++ (MPI/hybrid parallel)

.. code-block:: c

    pkomega_CG_R_update(v12, v2, x, r_l, status);
    pkomega_CG_C_update(v12, v2, x, r_l, status);
    pkomega_COCG_update(v12, v2, x, r_l, status);
    pkomega_BiCG_update(v12, v2, v14, v4, x, r_l, status);

Parameters

-  ``v12``

   ``DOUBLE PRECISION`` (for ``CG_R_update``), ``DOUBLE COMPLEX`` (for other cases). 
   The array with the length of ``ndim``. In/Output. 
   The product of the residual vector (``v2``) and the matrix.
   This routine returns the 2-norm of the updated residual vector
   as a first element of this array.
   This returned value is used, for examples, for printing the convergence profile.

-  ``v2``

   ``DOUBLE PRECISION`` (for ``CG_R_update``), ``DOUBLE COMPLEX`` (for other cases).
   The array with the length of ``ndim``. In/Output.
   The residual vector is input and the updated residual vector is output.

-  ``v14``

   Only used for ``BiCG_update``. ``DOUBLE COMPLEX``.
   The array with the length of ``ndim``. In/Output.
   The product of the shadow residual vector (``v4``) and the matrix is input.

-  ``v4``

   Only used for ``BiCG_update``. ``DOUBLE COMPLEX``.
   The array with the length of ``ndim``. In/Output.
   The shadow residual vector is input and the updated vector is output.

-  ``status``

   ``INTEGER``. The array with the length of ``3``. Output.
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
       ``?_init`` routine and the convergence condition is not satisfied
       at the ``itermax``\ -th iteration.
       ``2`` is returned if :math:`\alpha` diverges.
       ``3`` is returned if :math:`\pi_{\rm seed}` becomes 0.
       In the case of ``COCG_restart`` or ``BiCG_restart``,
       ``4`` is returned if the residual vector and
       the shadow residual vector are orthogonal.
       In other cases, ``0`` is returned.

   Third component (``status(3)``)
       The index of the seed point is returned.

``komega_????_getcoef``, ``pkomega_????_getcoef``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get the coefficients used in the restart calculation.
To call these routines,
``itermax`` in ``?_init`` routine must not be ``0`` .

Syntax

Fortran (Serial/OpenMP)

.. code-block:: fortran

    CALL komega_CG_R_getcoef(alpha_save, beta_save, z_seed, r_l_save)
    CALL komega_CG_C_getcoef(alpha_save, beta_save, z_seed, r_l_save)
    CALL komega_COCG_getcoef(alpha_save, beta_save, z_seed, r_l_save)
    CALL komega_BiCG_getcoef(alpha_save, beta_save, z_seed, r_l_save)

Fortran (MPI/hybrid parallel)

.. code-block:: fortran

    CALL pkomega_CG_R_getcoef(alpha_save, beta_save, z_seed, r_l_save)
    CALL pkomega_CG_C_getcoef(alpha_save, beta_save, z_seed, r_l_save)
    CALL pkomega_COCG_getcoef(alpha_save, beta_save, z_seed, r_l_save)
    CALL pkomega_BiCG_getcoef(alpha_save, beta_save, z_seed, r_l_save)

C/C++ (Serial/OpenMP)

.. code-block:: c

    komega_CG_R_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
    komega_CG_C_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
    komega_COCG_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
    komega_BiCG_getcoef(alpha_save, beta_save, &z_seed, r_l_save);

C/C++ (MPI/hybrid parallel)

.. code-block:: c

    pkomega_CG_R_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
    pkomega_CG_C_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
    pkomega_COCG_getcoef(alpha_save, beta_save, &z_seed, r_l_save);
    pkomega_BiCG_getcoef(alpha_save, beta_save, &z_seed, r_l_save);

Parameters

-  ``alpha_save``

   ``DOUBLE PRECISION`` (for ``CG_R_getoef``, ``CG_C_getoef``),
   ``DOUBLE COMPLEX`` (for other cases).
   The array with the length of the number of maximum iteration. Output.
   The parameters :math:`\alpha` of the (Bi)CG method at each iteration.

-  ``beta_save``

   ``DOUBLE PRECISION`` (for ``CG_R_getoef``, ``CG_C_getoef``),
   ``DOUBLE COMPLEX`` (for other cases).
   The array with the length of the number of maximum iteration. Output.
   The parameters :math:`\beta` of the (Bi)CG method at each iteration.

-  ``z_seed``

   ``DOUBLE PRECISION`` (for ``CG_R_getoef``, ``CG_C_getoef``),
   ``DOUBLE COMPLEX`` (for other cases). Scalar. Output. Seed shift.

-  ``r_l_save``

   ``DOUBLE PRECISION`` (for ``CG_R_getoef``), ``DOUBLE COMPLEX`` (for other cases).
   The array with the length of the number of maximum iteration :math:`\times` ``nl``.
   Output.
   The projected residual vectors at each iteration.

``komega_????_getvec``, ``pkomega_????_getvec``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get the residual vectors to use the restart calculation.
To call these routines,
``itermax`` in the ``?_init`` routine must not be ``0``.

Syntax

Fortran (Serial/OpenMP)

.. code-block:: fortran

    CALL komega_CG_R_getvec(r_old)
    CALL komega_CG_C_getvec(r_old)
    CALL komega_COCG_getvec(r_old)
    CALL komega_BiCG_getvec(r_old, r_tilde_old)

Fortran (MPI/hybrid parallel)

.. code-block:: fortran

    CALL pkomega_CG_R_getvec(r_old)
    CALL pkomega_CG_C_getvec(r_old)
    CALL pkomega_COCG_getvec(r_old)
    CALL pkomega_BiCG_getvec(r_old, r_tilde_old)

C/C++ (Serial/OpenMP)

.. code-block:: c

    komega_CG_R_getvec(r_old);
    komega_CG_C_getvec(r_old);
    komega_COCG_getvec(r_old);
    komega_BiCG_getvec(r_old, r_tilde_old);

C/C++ (MPI/hybrid parallel)

.. code-block:: c

    pkomega_CG_R_getvec(r_old);
    pkomega_CG_C_getvec(r_old);
    pkomega_COCG_getvec(r_old);
    pkomega_BiCG_getvec(r_old, r_tilde_old);

Parameters

-  ``r_old``

   ``DOUBLE PRECISION`` (for ``CG_R_getvec``), ``DOUBLE COMPLEX`` (for other cases).
   The array with the length of ``ndim``. Output.
   The residual vector at the second last step in the previous calculation.

-  ``r_tilde_old``

   Only used for ``BiCG_getvec``. ``DOUBLE COMPLEX``. 
   The array with the length of ``ndim``. Output.
   The shadow residual vector at the second last step in the previous calculation.

``komega_????_getresidual``, ``pkomega_????_getresidual``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Get the values of 2-norm of the residual vector at each shift points.
These routines can be called from anywhere between ``komega_????_init``
and ``komega_????_finalize``.
These routines do not affect the calculation results.

Syntax

Fortran (Serial/OpenMP)

.. code-block:: fortran

    CALL komega_CG_R_getresidual(res)
    CALL komega_CG_C_getresidual(res)
    CALL komega_COCG_getresidual(res)
    CALL komega_BiCG_getresidual(res)

Fortran (MPI/hybrid parallel)

.. code-block:: fortran

    CALL pkomega_CG_R_getresidual(res)
    CALL pkomega_CG_C_getresidual(res)
    CALL pkomega_COCG_getresidual(res)
    CALL pkomega_BiCG_getresidual(res)

C/C++ (Serial/OpenMP)

.. code-block:: c

    komega_CG_R_getresidual(res);
    komega_CG_C_getresidual(res);
    komega_COCG_getresidual(res);
    komega_BiCG_getresidual(res);

C/C++ (MPI/hybrid parallel)

.. code-block:: c

    pkomega_CG_R_getresidual(res);
    pkomega_CG_C_getresidual(res);
    pkomega_COCG_getresidual(res);
    pkomega_BiCG_getresidual(res);

Parameters

-  ``res``

   ``DOUBLE PRECISION``. The array with the length of ``nz``. Output.
   The values of 2-norm of the residual vector at each shift points are
   returned.

``komega_????_finalize``, ``pkomega_????_finalize``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Release memories of the arrays stored in the library.

Syntax

Fortran (Serial/OpenMP)

.. code-block:: fortran

    CALL komega_CG_R_finalize()
    CALL komega_CG_C_finalize()
    CALL komega_COCG_finalize()
    CALL komega_BiCG_finalize()

Fortran (MPI/hybrid parallel)

.. code-block:: fortran

    CALL pkomega_CG_R_finalize()
    CALL pkomega_CG_C_finalize()
    CALL pkomega_COCG_finalize()
    CALL pkomega_BiCG_finalize()

C/C++ (Serial/OpenMP)

.. code-block:: c

    komega_CG_R_finalize();
    komega_CG_C_finalize();
    komega_COCG_finalize();
    komega_BiCG_finalize();

C/C++ (MPI/hybrid parallel)

.. code-block:: c

    pkomega_CG_R_finalize();
    pkomega_CG_C_finalize();
    pkomega_COCG_finalize();
    pkomega_BiCG_finalize();

Sample codes for using shifted BiCG library
-------------------------------------------

As a typical example, the usage of shifted BiCG library is shown below.

.. code-block:: fortran

    PROGRAM my_prog
      !
      USE komega_bicg, ONLY : komega_BiCG_init, komega_BiCG_restart, &
      &                       komega_BiCG_update, komega_BiCG_getcoef, &
      &                       komega_BiCG_getvec, komega_BiCG_finalize
      USE solve_cc_routines, ONLY : input_size, input_restart, &
      &                             projection, &
      &                             hamiltonian_prod, generate_system, &
      &                             output_restart, output_result
      !
      IMPLICIT NONE
      !
      INTEGER,SAVE :: &
      & rnd_seed, &
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
      & itermin, & ! First iteration in this run
      & iter,    & ! Counter for Iteration
      & status(3)
      !
      ! Input Size of vectors
      !
      CALL input_size(ndim,nl,nz)
      !
      ALLOCATE(v12(ndim), v2(ndim), v14(ndim), v4(ndim), r_l(nl), &
      &        x(nl,nz), z(nz), ham(ndim,ndim), rhs(ndim))
      !
      CALL generate_system(ndim, ham, rhs, z)
      !
      ! Check: Whether the restart file is exist.
      !
      CALL input_restart(iter_old, zseed, alpha, beta, r_l_save)
      !
      WRITE(*,*)
      WRITE(*,*) "#####  CG Initialization  #####"
      WRITE(*,*)
      !
      IF(iter_old > 0) THEN
        !
        ! When restarting, counter
        !
        itermin = iter_old + 1
        CALL komega_BiCG_restart(ndim, nl, nz, x, z, max(0,itermax), &
        &                        threshold, &
        &                 status, iter_old, v2, v12, v4, v14, alpha, &
        &                 beta, z_seed, r_l_save)
        !
        ! These vectors were saved in BiCG routine
        !
        DEALLOCATE(alpha, beta, r_l_save)
        !
        IF(status(1) /= 0) GOTO 10
        !
      ELSE
         !
         itermin = 1
         !
         ! Generate Right Hand Side Vector
         !
         v2(1:ndim) = rhs(1:ndim)
         v4(1:ndim) = CONJG(v2(1:ndim))
         !v4(1:ndim) = v2(1:ndim)
         !
         CALL komega_BiCG_init(ndim, nl, nz, x, z, max(0,itermax), &
         &                     threshold)
         !
      END IF
      !
      ! BiCG Loop
      !
      WRITE(*,*)
      WRITE(*,*) "#####  CG Iteration  #####"
      WRITE(*,*)
      !
      DO iter = 1, abs(itermax)
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
         CALL komega_BiCG_update(v12, v2, v14, v4, x, r_l, status)
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
      iter_old = ABS(status(1))
      !
      ! Get these vectors for restart in the Next run
      !
      IF(itermax > 0) THEN
         !
         ALLOCATE(alpha(iter_old), beta(iter_old), r_l_save(nl,iter_old))
         !
         CALL komega_BiCG_getcoef(alpha, beta, z_seed, r_l_save)
         CALL komega_BiCG_getvec(v12,v14)
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
      CALL komega_BiCG_finalize()
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

