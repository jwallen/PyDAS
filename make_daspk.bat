@echo off

:: Set the Fortran compiler and compiler flags
set f77=gfortran
set cflags=-O3

echo Compiling DASPK3.1...
CALL %f77% %cflags% -c daspk31/solver/adf_dummy.f -o daspk/solver/adf_dummy.o
CALL %f77% %cflags% -c daspk31/solver/daux.f -o daspk/solver/daux.o
CALL %f77% %cflags% -c daspk31/solver/ddaskr.f -o daspk/solver/ddaskr.o
CALL %f77% %cflags% -c daspk31/solver/ddaspk.f -o daspk/solver/ddaspk.o
CALL %f77% %cflags% -c daspk31/solver/dlinpk.f -o daspk/solver/dlinpk.o
CALL %f77% %cflags% -c daspk31/solver/dsensd.f -o daspk/solver/dsensd.o
CALL %f77% %cflags% -c daspk31/solver/mpi_dummy.f -o daspk/solver/mpi_dummy.o

CALL ar rcs dassl/libddaspk31.a daspk31/solver/adf_dummy.o daspk31/solver/daux.o daspk31/solver/ddaskr.o daspk31/solver/ddaspk.o daspk31/solver/dlinpk.o daspk31/solver/dsensd.o daspk31/solver/mpi_dummy.o

echo Compiling PyDAS DASPK wrapper...
CALL python setup.py build_ext daspk --compiler=mingw32 --inplace

:end
pause
