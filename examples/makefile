examples ?=
root ?= ../

include $(root)makefile.in

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

bvp.o: $(examples)bvp.f90
	$(FC) -o $(examples)bvp.o -c $(examples)bvp.f90 -J$(examples)../include -I/usr/include

bvp: $(examples)bvp.o $(LIB)liblinalg.a $(LIB)libmesh.a $(LIB)libfem.a
	$(FC) -o $(examples)bvp $(examples)bvp.o -llinalg -lmesh -lfem -llapack -lblas -lnetcdf -lnetcdff -L$(LIB)

f2k3_to_cxx.o: $(examples)f2k3_to_cxx.cpp
	$(CXX) -c $(examples)f2k3_to_cxx.cpp -I$(INCLUDE) -I/home/daniel/programs/tetgen1.4.3

f2k3_to_cxx: $(examples)f2k3_to_cxx.o
	$(CXX) -o $(examples)f2k3_to_cxx $(examples)f2k3_to_cxx.o -llinalg -llapack -lblas -ltet -L$(LIB) -L/home/daniel/programs/tetgen1.4.3

wiffle.1.node: $(examples)wiffle.mesh
	tetgen -pne wiffle.mesh

