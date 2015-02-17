################################################################################
#
#   Makefile for DASPK3.1
#
################################################################################

F77=gfortran

CFLAGS=-fPIC -O3

DOBJ=solver/adf_dummy.o solver/daux.o solver/ddaskr.o solver/ddaspk.o solver/dlinpk.o solver/dsensd.o solver/mpi_dummy.o

DLIB=libddaspk31.a

all: $(DLIB)

$(DLIB): $(DOBJ)
	ar rcs $(DLIB) $(DOBJ)

%.o: %.f
	$(F77) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(DLIB) $(DOBJ)
