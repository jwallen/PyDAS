*************************************************************************
PyDAS - A Python wrapper to several differential algebraic system solvers
*************************************************************************

Introduction
============

PyDAS provides a means for Python code to utilize several notable Fortran-based
differential algebraic system solvers from Python code. The solvers made
available -- DASSL, DASPK, and DASKR -- are all publicly-available from 
`Netlib <http://www.netlib.org/ode/>`_, and are distributed with PyDAS. PyDAS
provides a Python extension type for each solver, which in turn provides a
Pythonic means of setting the solver options, providing residual and jacobian
functions, and running the solver.

The DASSL, DASPK, and DASKR solvers are all substantially more robust than
VODE, the solver used within the ODE solver functionality provided by 
`SciPy <http://www.scipy.org/>`_.

License
=======

Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu).

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the 'Software'),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

Dependencies
============

PyDAS has been tested on Python versions 2.5 and 2.6. It may or may not work
for other Python versions.

There are several Python dependencies that you must install before installing 
PyDAS:

* Cython (0.12 or later)
* NumPy (1.3 or later)

In addition, you will also need a Fortran compiler and a C compiler that
produce object files that can interoperate. The ``gfortran`` and ``gcc`` 
compiles from the GNU Compiler Collection are known to work.

The code for the differential algebraic system solvers DASSL, DASPK, and DASKR
has been provided with the PyDAS package. The licenses for these solvers is
different than that of the PyDAS wrapper code. **You are responsible for knowing
and abiding by all licenses associated with each solver as well as with PyDAS
as a whole.**

Installation
============

A Makefile has been provided that can be used to compile all of the solvers
and the PyDAS wrapper code. To use, invoke the following command from the
base package directory::

	$ make

You may wish to write a file `make.inc` that sets certain variables used by
the Makefiles (e.g. the Fortran compiler).

