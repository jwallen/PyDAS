#!/usr/bin/python
# -*- coding: utf-8 -*-

if __name__ == '__main__':
	
	from distutils.core import setup
	from distutils.extension import Extension
	from Cython.Distutils import build_ext
	
	import Cython.Compiler
	Cython.Compiler.Options.annotate = True
	
	# The Cython modules to setup
	ext_modules = [
		Extension('dassl', ['dassl.pyx'], library_dirs=['dassl'], libraries=['gfortran', 'ddassl']),
	]

	setup(cmdclass = {'build_ext': build_ext},
		ext_modules = ext_modules
	)