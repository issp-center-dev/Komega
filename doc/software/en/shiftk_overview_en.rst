Overview
========

This is document for the sample program which uses
K :math:`\omega` library in the ISSP Math Library;
this program computes the Green's function with :math:`K(\omega)`.
For the details of K :math:`\omega` library, See
":math:`K(\omega)` manual" in this package.

Calculation in this program
---------------------------

This program compute the Green's function

.. math::

   \begin{align}
   G_{i}(z) =
   \langle i | (z-{\hat H})^{-1}| i \rangle
   \equiv 
   {\boldsymbol \varphi}_i^{*} \cdot (z-{\hat H})^{-1} {\boldsymbol \varphi}_i,
   \end{align}

where :math:`| i \rangle` is a wavefunction,
:math:`{\cal H}` is the Hamiltonian, and
:math:`z` is a complex frequency.

:math:`{\cal H}` in the above equation is obtained by either the
      following two ways:

-  Input :math:`{\cal H}` as a file with the MatrixMarket format

-  Construct :math:`{\cal H}` as a Hamiltonian of the
   Heisenberg model in this program.

In the computation of the Green's function,
we use either the following two method according to the type
of :math:`{\hat H}` (a real- or a complex- number).

-  :math:`{\hat H}` of real numbers : Shifted Bi-Conjugate Gradient(BiCG) method

-  :math:`{\hat H}` of complex numbers : Shifted Conjugate Orthogonal Conjugate Gradient(COCG) method

