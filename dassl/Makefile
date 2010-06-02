################################################################################
#
#   Makefile for DASSL
#
################################################################################

F77=gfortran

CFLAGS=-fPIC

OBJ=daux.o ddassl.o dlinpk.o

LIB=libddassl.a

all: dassl

dassl: $(OBJ)
	ar rcs $(LIB) $(OBJ)

%.o: %.f
	$(F77) $(CFLAGS) -c $< -o $@
	
clean:
	rm $(LIB) $(OBJ)
