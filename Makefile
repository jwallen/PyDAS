################################################################################
#
#   Makefile for PyDAS
#
################################################################################

CYTHON_FLAGS=--inplace

CLEAN_FLAGS=

-include make.inc

.PHONY: all build install clean

all: build

build:
	python setup.py build_ext $(CYTHON_FLAGS)

install:
	python setup.py install

clean:
	python setup.py clean $(CLEAN_FLAGS)
	rm -f *.so *.pyc dassl/_dassl-f2pywrappers.f dassl/_dasslmodule.c daspk/_daspk-f2pywrappers.f daspk/_daspkmodule.c daskr/_daskr-f2pywrappers.f daskr/_daskrmodule.c
