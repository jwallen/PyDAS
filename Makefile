################################################################################
#
#   Makefile for PyDAS
#
################################################################################

-include make.inc

all: DASSL DASPK DASKR cython

cython:
	python setup.py build_ext $(CYTHON_FLAGS)

install:
	python setup.py install

DASSL:
	$(MAKE) -C dassl F77=$(F77)

DASPK:
	$(MAKE) -C daspk F77=$(F77)

DASKR:
	$(MAKE) -C daskr F77=$(F77)

clean: clean-DASSL clean-DASPK clean-DASKR clean-cython

clean-DASSL:
	$(MAKE) -C dassl clean

clean-DASPK:
	$(MAKE) -C daspk clean

clean-DASKR:
	$(MAKE) -C daskr clean

clean-cython:
	python setup.py clean $(CLEAN_FLAGS)
	rm -f *.so *.pyc

help:
	@echo ""
	@echo "This makefile can be used to build PyDAS and its dependencies."
	@echo ""
	@echo "Typing \`make\` with no arguments will compile all three DAE solvers (DASSL,"
	@echo "DASPK, and DASKR) to static libraries and compile the PyDAS Python modules"
	@echo "that provide the Python interface to these solvers."
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
	@echo "    DASKR    to compile the DASKR solver"
	@echo "    cython   to compile the PyDAS Python wrapper modules"
	@echo ""

