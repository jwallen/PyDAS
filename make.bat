@echo off

:: Set the Fortran compiler and compiler flags
set f77=gfortran
set cflags=-O3

echo Compiling DASSL...
CALL %f77% %cflags% -c dassl/daux.f -o dassl/daux.o
CALL %f77% %cflags% -c dassl/ddassl.f -o dassl/ddassl.o
CALL %f77% %cflags% -c dassl/dlinpk.f -o dassl/dlinpk.o
CALL ar rcs dassl/libddassl.a dassl/daux.o dassl/ddassl.o dassl/dlinpk.o

echo Compiling PyDAS...
CALL python setup.py build_ext --compiler=mingw32 --inplace

:end
pause
