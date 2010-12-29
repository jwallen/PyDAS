@echo off

:: Set the Fortran compiler and compiler flags
set f77=__gfortran
set cflags=-O3

echo Compiling DASSL...
%f77% %cflags% -c dassl/daux.f -o dassl/daux.o
%f77% %cflags% -c dassl/ddassl.f -o dassl/ddassl.o
%f77% %cflags% -c dassl/dlinpk.f -o dassl/dlinpk.o
ar rcs dassl/libddassl.a dassl/daux.o dassl/ddassl.o dassl/dlinpk.o

echo Compiling PyDAS...
python setup.py build_ext --compiler=mingw32 --inplace

:end
pause
