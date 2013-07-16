examples ?=
include $(examples)../makefile.in

INCLUDE ?= $(examples)../include/
LIB ?= $(examples)../lib/

bvp.o: $(examples)bvp.f90
	$(FC) -o $(examples)bvp.o -c $(examples)bvp.f90 -J$(examples)../include -I/usr/include

bvp: $(examples)bvp.o $(LIB)liblinalg.a $(LIB)libmesh.a $(LIB)libfem.a
	$(FC) -o $(examples)bvp $(examples)bvp.o -llinalg -lmesh -lfem -llapack -lblas -lnetcdf -lnetcdff -L$(LIB)

