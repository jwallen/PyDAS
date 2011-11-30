#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   PyDAS - A Python wrapper to several differential algebraic system solvers
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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

try:
    import numpy
except ImportError:
    print('The numpy package is required to install PyDAS.')
    quit()

from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup

################################################################################

def configuration(parent_package='',top_path=None):
    """
    Return information about the extension modules for f2py to build.
    """
    config = Configuration('pydas', parent_package, top_path)
    config.add_extension('_dassl', ['dassl/ddassl.f', 'dassl/daux.f', 'dassl/dlinpk.f', 'dassl/ddassl.pyf'])
    config.add_extension('_daspk', ['daspk/solver/daux.f', 'daspk/solver/ddaspk.f', 'daspk/solver/dlinpk.f', 'daspk/preconds/dbanpre.f', 'daspk/preconds/dilupre.f', 'daspk/preconds/drbdpre.f', 'daspk/preconds/drbgpre.f', 'daspk/preconds/dsparsk.f', 'daspk/ddaspk.pyf'])
    config.add_extension('_daskr', ['daskr/solver/daux.f', 'daskr/solver/ddaskr.f', 'daskr/solver/dlinpk.f', 'daskr/preconds/dbanpre.f', 'daskr/preconds/dilupre.f', 'daskr/preconds/drbdpre.f', 'daskr/preconds/drbgpre.f', 'daskr/preconds/dsparsk.f', 'daskr/ddaskr.pyf'])
    return config

setup(
    name = 'PyDAS',
    version = '0.1.0',
    description = 'A Python wrapper to several differential algebraic system solvers',
    author = 'Joshua W. Allen',
    author_email = 'jwallen@mit.edu',
    url = 'http://github.com/jwallen/PyDAS',
    packages = ['pydas'],
    configuration = configuration,
    requires = ['numpy (>=1.3.0)'],
    provides = ['pydas'],
)
