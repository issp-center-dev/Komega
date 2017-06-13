Re-distribution of this library
===============================

Contain Kω in your program
---------------------------

:math:`K\omega` library is distributed with the :ref:`lgplicense` (LGPL).
It is summarized as follows:

 * :math:`K\omega` can be freely distributed, modified, copied and pasted,
   in a private program (in the research group, co-workers, etc.).
   
 * For the released program (open-source, free, commercial software etc.):
   
    * When you contain the source-code of :math:`K\omega` (either as is and modified)
      in the distributed source code of your program,
      please distribute your program with LGPL/GPL.
      
    * If you do not include the source-code of :math:`K\omega` (just call it),
      you can freely distribute your program with any licenses.
      
    * If you distribute a binary file which is statically linked to :math:`K\omega` library,
      please use LGPL/GPL. However, if you distribute a binary file which is dynamically linked to
      :math:`K\omega` library (therefore :math:`K\omega` itself is not contained),
      you can freely distribute your binary file with any licenses.

Build Kω without Autoconf
--------------------------

In this package, :math:`K\omega` is built with Autotools (Autoconf, Automake, Libtool).
If you do not want to use Autotools for your distributed program with :math:`K\omega` source,
you can use the following simple Makefile (please care about TAB).

.. code-block:: makefile

   F90 = gfortran
   FFLAGS = -fopenmp -g -O2 #-D__MPI -D__NO_ZDOT -D__KOMEGA_THREAD
   
   .SUFFIXES :
   .SUFFIXES : .o .F90
   
   OBJS = \
   komega_cg_c.o \
   komega_cg_r.o \
   komega_cocg.o \
   komega_bicg.o \
   komega_math.o \
   komega_vals.o
   
   all:libkomega.a
   
   libkomega.a:$(OBJS)
        ar cr libkomega.a $(OBJS)
   
   .F90.o:
        $(F90) -c $< $(FFLAGS)
   
   clean:
        rm -f *.o *.a *.mod
   
   komega_cg_c.o:komega_math.o
   komega_cg_c.o:komega_vals.o
   komega_cg_r.o:komega_math.o
   komega_cg_r.o:komega_vals.o
   komega_cocg.o:komega_math.o
   komega_cocg.o:komega_vals.o
   komega_bicg.o:komega_math.o
   komega_bicg.o:komega_vals.o
   komega_math.o:komega_vals.o

Preprocessor macros ``__MPI``, ``__NO_ZDOT``, and ``__KOMEGA_THREAD`` correspond to
``--with-mpi=yes``, ``--disable-zdot``, and ``--enable-thread`` of the options of ``configure``, respectively.
      
.. _lgplicense:
      
Lesser General Public License
-----------------------------

*© 2016- The University of Tokyo. All rights reserved.*

| This software is developed under the support of
| "*Project for advancement of software usability in materials science*" by The
| Institute for Solid State Physics, The University of Tokyo.
|
| This library is free software; you can redistribute it and/or
| modify it under the terms of the GNU Lesser General Public
| License as published by the Free Software Foundation; either
| version 2.1 of the License, or (at your option) any later version.
| This library is distributed in the hope that it will be useful,
| but WITHOUT ANY WARRANTY; without even the implied warranty of
| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
| Lesser General Public License for more details.
|
| You should have received a copy of the GNU Lesser General Public
| License along with this library; if not, write to the Free Software
| Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
|
| For more details, See 'COPYING.LESSER' in the root directory of this library.

