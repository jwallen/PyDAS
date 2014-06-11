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

* `Python <http://www.python.org/>`_ (versions 2.5.x and 2.6.x are known to work)

* `NumPy <http://numpy.scipy.org/>`_ (version 1.3.0 or later is recommended)

* `Cython <http://www.cython.org/>`_ (version 0.12.1 or later is recommended)

In addition, you will also need a Fortran compiler and a C compiler that
produce object files that can interoperate. The ``gfortran`` and ``gcc`` 
compiles from the GNU Compiler Collection are known to work. On Windows the
`MinGW <http://www.mingw.org/>`_ compiler collection provides these compilers.

The code for the differential algebraic system solvers DASSL, DASPK, and DASKR
has been provided with the PyDAS package. The licenses for these solvers is
different than that of the PyDAS wrapper code. **You are responsible for knowing
and abiding by all licenses associated with each solver as well as with PyDAS
as a whole.**

Installation
============

.. note:: 

    Currently only the DASSL solver has been wrapped. The installation 
    scripts therefore only build and install the DASSL wrapper by default.

Windows
-------

The provided batch scripts will compile all of the solvers and the PyDAS
wrapper code. These scripts presume that you have the 32-bit version of the
MinGW C and Fortran compilers installed. Once you have run the batch script,
you can install PyDAS into your Python packages if you desire by running the
following command from the base package directory:

    > python setup.py install

Linux
-----

A Makefile has been provided that can be used to compile all of the solvers
and the PyDAS wrapper code. To use, invoke the following command from the
base package directory::

    $ make

This command will build PyDAS in-place, rather than installing it to your
Python package directory. If you wish to formall install PyDAS, run the
following command from the base package directory after the ``make`` command
(you may need root privileges for this)::

    $ python setup.py install

You may wish to write a file `make.inc` that sets certain variables used by
the Makefiles (e.g. the Fortran compiler). An example of such a file, 
`make.inc.example`, has been provided.


Mac OS X
---------

Homebrew (http://brew.sh) is an easy way to get gfortran::

    $ brew install gcc

But your system may still not be able to find the correct `libgfortran.a` library file.
This should make it work::

    $ export LIBRARY_PATH=$(dirname `gfortran -print-libgcc-file-name`)
    $ make F77=gfortran

Then, to get it installed into your proper python place (you may need to prefix this with ``sudo ``)::

    $ make install F77=gfortran
