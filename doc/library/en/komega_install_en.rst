Install
=======

Overall procedure
-----------------

First, please type

.. code-block:: bash

   $ ./configure --prefix=install_dir

Then, this script checks the compiler and the libraries required for the installation,
and creates Makefiles.
``install_dir`` indicates the full path of the directory where the library is installed
(you should replace it according to your case).
If none is specified, ``/use/local/`` is chosen for storing libraries
by ``make install``  (Therefore, if one is not the admin, ``install_dir`` must be specified to
the different directory).
``configure`` has many options, and they are used according to the environment etc.
For more details, please see :ref:`configoption`.

After ``configure`` finishes successfully and Makefiles are generated,
please type

.. code-block:: bash

   $ make

to build libraries.

You can optionally run the regression tests with

.. code-block:: bash

   $ make check

which runs the sample solvers and checks that they converge to within a
tolerance (it returns a non-zero status on failure).

Then please type

.. code-block:: bash

   $ make install

to store libraries and the sample program to ``install_dir/lib`` and ``install_dir/bin``, respectively.
Although one can use libraries and the sample program without ``make install``,
they are a little different to the installed one.

Add the :math:`K\omega` library directory (``install_dir/lib``) to the
search path of the dynamically linked program (environment variable ``LD_LIBRARY_PATH``).

.. code-block:: bash

   $ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:install_dir/lib

.. _configoption:

Options for configure
---------------------

``configure`` has many options and environment variables.
They can be specified at once. E.g.

.. code-block:: bash

  $ ./configure --prefix=/home/komega/ --with-mpi=yes FC=mpif90

All options and variables have default values.
We show a part of them as follows:

``---prefix``

   Default: ``---prefix=/usr/local/``.
   Specify the directory where the library etc. are installed.

``--with-mpi``

   Default: ``--with-mpi=no`` (without MPI).
   Whether use MPI (``--with-mpi=yes``), or not.

``--with-openmp``

   Default: ``--with-openmp=yes`` (with OpenMP).
   Whether use OpenMP or not (``--with-openmp=no``).

``--enable-shared``

   Default: ``--enable-shared``.
   Whether generate shared library.

``--enable-static``

   Default: ``--enable-static``.
   Whether generate static library.

``--enable-zdot`` / ``--disable-zdot``

   Default: autodetection.
   Controls whether the BLAS ``ZDOTC`` / ``ZDOTU`` functions are used.
   With some BLAS libraries (notably the standard BLAS and the Accelerate
   framework on macOS) these functions return their complex result through a
   hidden first argument (the ``f2c``/``g77`` calling convention) instead of in
   registers as :math:`K\omega` expects, which makes the program crash.
   By default ``configure`` autodetects this by compiling and running a small
   test program: if the BLAS ``ZDOTC`` works it is used, otherwise the library
   falls back to the Fortran intrinsics ``DOT_PRODUCT`` / ``SUM`` (equivalent to
   defining the ``-D__NO_ZDOT`` macro).
   Pass ``--enable-zdot`` to force the use of the BLAS functions, or
   ``--disable-zdot`` to force the intrinsic fallback; either choice bypasses
   the autodetection.
   When cross-compiling (the test program cannot be run) the safe fallback
   (``--disable-zdot``) is chosen automatically.

``--enable-threadsafe``

   Default: ``--disable-threadsafe``.
   If you want to call :math:`K\omega` routine in the parallel region
   (i.e. plan to solve different equations among threads),
   please use this option (**Experimental**).

``FC``

   Default: The fortran compiler chosen automatically from those in the system.
   When ``--with-mpi`` is specified, the corresponding MPI compiler
   (such as ``mpif90``) is searched.
   If ``FC`` printed the end of the standard-output of ``configure`` is not
   what you want, please set it manually as ``./configure FC=gfortran``.

``--help``

   Display all options including above, and stop without configuration.
