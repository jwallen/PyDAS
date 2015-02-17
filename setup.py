#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   PyDAS - A Python wrapper to several differential algebraic system solvers
#
#   Copyright (c) 2010-2014 by Joshua W. Allen (joshua.w.allen@gmail.com) 
#                           extended by Connie W. Gao (connieg@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import numpy
import sys

if __name__ == '__main__':
    
    # Use setuptools by default (requires Python>=2.6) if available
    # If not available, fall back to distutils
    # Using setuptools enables support for compiling wheels
    try:
        from setuptools import setup, Extension
    except ImportError:
        from distutils.core import setup
        from distutils.extension import Extension

    from Cython.Distutils import build_ext
    
    # Turn on HTML annotation file generation
    import Cython.Compiler.Options
    Cython.Compiler.Options.annotate = True
    
    # The Cython extension modules to compile
    pydas_ext = Extension(
            'pydas.dassl', 
            ['pydas/dassl.pyx'], 
            include_dirs=['pydas', numpy.get_include()], 
            libraries=['gfortran'], 
            extra_objects=['dassl/daux.o','dassl/ddassl.o','dassl/dlinpk.o'],
        )
    pydaspk_ext = Extension(
            'pydas.daspk', 
            ['pydas/daspk.pyx'], 
            include_dirs=['pydas', numpy.get_include()], 
            libraries=['gfortran'], 
            extra_objects=['daspk31/solver/adf_dummy.o','daspk31/solver/daux.o','daspk31/solver/ddaspk.o','daspk31/solver/dlinpk.o','daspk31/solver/dsensd.o','daspk31/solver/mpi_dummy.o'],
        )
    

    modules = ['pydas.dassl']
    extensions = [pydas_ext]

    if 'daspk' in sys.argv:
        # Optionally compile and make pydaspk if the user requests it
        sys.argv.remove('daspk')
        modules.append('pydas.daspk')
        extensions.append(pydaspk_ext)

    # Run the setup command
    setup(name='PyDAS',
        version='0.1.0',
        description='A Python wrapper to several differential algebraic system solvers',
        author='Joshua W. Allen',
        author_email='joshua.w.allen@gmail.com',
        url='http://github.com/jwallen/PyDAS',
        py_modules= modules,
        packages = ['pydas'],
        package_data = {'pydas': ['*.pxd']},
        cmdclass = {'build_ext': build_ext},
        ext_modules = extensions
    )
