Overview
========

This document is a manual for :math:`K(\omega)` which is the library to
solve the shifted linear equation within the Krylov subspace.
This library provides routines to solve the following shifted linear equation
(with the projection),

.. _shiftedeq:

.. math::

   \begin{align}
     G_{i j}(z) = \langle i | (z {\hat I} -{\hat H})^{-1}| j \rangle \equiv 
     {\boldsymbol \varphi}_i^{*} \cdot (z{\hat I}-{\hat H})^{-1} {\boldsymbol \varphi}_j.
     \end{align}

The source codes of :math:`K(\omega)` is written in FORTRAN
and requires the BLAS Level 1 routines.

