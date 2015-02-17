################################################################################
#
#   Makefile for PyDAS
#
################################################################################

F77=gfortran

CYTHON_FLAGS=--inplace

-include make.inc

.PHONY: DASSL DASPK DASPK31 DASKR cython clean 

all: DASSL DASPK DASKR cython

daspk: DASPK31 cython-daspk

cython-daspk:
	python setup.py build_ext daspk $(CYTHON_FLAGS)

cython: DASSL DASPK DASKR 
	python setup.py build_ext $(CYTHON_FLAGS)

install: DASSL DASPK DASKR cython
	python setup.py install

DASSL:
	$(MAKE) -C dassl F77=$(F77)

DASPK:
	$(MAKE) -C daspk F77=$(F77)

DASPK31:
	$(MAKE) -C daspk31 F77=$(F77)

DASKR:
	$(MAKE) -C daskr F77=$(F77)

clean: clean-DASSL clean-DASPK clean-DASPK31 clean-DASKR clean-cython
	rm -rf build

clean-DASSL:
	$(MAKE) -C dassl clean

clean-DASPK:
	$(MAKE) -C daspk clean

clean-DASPK31:
	$(MAKE) -C daspk31 clean

clean-DASKR:
	$(MAKE) -C daskr clean

clean-cython:
	python setup.py clean $(CLEAN_FLAGS)
	rm -f pydas/*.so pydas/*.pyc pydas/*.c

help:
	@echo ""
	@echo "This makefile can be used to build PyDAS and its dependencies."
	@echo ""
	@echo "Typing \`make\` with no arguments will compile all three DAE solvers (DASSL,"
	@echo "DASPK, and DASKR) to static libraries and compile the PyDAS Python modules"
	@echo "that provide the Python interface to these solvers."
	@echo ""
	@echo "Typing \`make daspk\` after typing \`make\` will then additionally compile 
	@echo "the optional DASPK 3.1 solver as well as the cython module pydaspk associated with it."
	@echo "The DASPK 3.1 fortran source files must first be downloaded externally and placed"
	@echo "in the daspk31 folder."
	@echo ""
	@echo "Typing \`make clean\` will delete all of the intermediate build files,"
	@echo "compiled libraries, and compiled Python modules for all three DAE solvers and"
	@echo "the PyDAS modules."
	@echo ""
	@echo "Individual dependencies can be specified using \`make <target>\`, where"
	@echo "<target> is one of:"
	@echo ""
	@echo "    DASSL    to compile the DASSL solver"
	@echo "    DASPK    to compile the DASPK solver"   
	@echo "    DASPK31    to compile the DASPK31 solver"   
	@echo "    DASKR    to compile the DASKR solver"
	@echo "    cython   to compile the PyDAS Python wrapper module for DASSL"
	@echo "    cython-daspk   to compile the PyDAS Python wrapper module for both DASSL and DASPK3.1"
	@echo ""

