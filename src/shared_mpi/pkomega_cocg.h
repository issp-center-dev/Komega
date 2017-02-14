/*
ISSP Math Library - A library for solving linear systems in materials science
Copyright (C) 2016 Mitsuaki Kawamura

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

For more details, See ÅeCOPYING.LESSERÅf in the root directory of this library.
*/

#include <mpi.h>
#include <complex.h>
#pragma once

void pkomega_cocg_init(int *ndim, int *nl, int *nz, double complex *x,
                       double complex *z, int *itermax, double *threshold,
                       MPI_Comm *comm);
void pkomega_cocg_restart(int *ndim, int *nl, int *nz, double complex *x,
                          double complex *z, int *itermax, double *threshold,
                          MPI_Comm *comm, int *status, int *iter_old, double complex *v2,
                          double complex *v12, 
                          double complex *alpha_save, double complex *beta_save,
                          double complex *z_seed, double complex *r_l_save);
void pkomega_cocg_update(double complex *v12, double complex *v2, double complex *x,
                         double complex *r_l, int *status);
void pkomega_cocg_getcoef(double complex *alpha_save, double complex *beta_save,
                          double complex *z_seed, double complex *r_l_save);
void pkomega_cocg_getvec(double complex *r_old);
void pkomega_cocg_getresidual(double *res);
void pkomega_cocg_finalize();
