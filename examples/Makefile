OBJECTS = ../lib/liblinalg.a ../lib/libmesh.a ../lib/libfem.a
LFLAGS = -lnetcdf -lnetcdff
IFLAGS = -I/usr/local/include -I../include
DFLAGS = -g -fbounds-check -Wall

FC = gfortran

.PHONY: clean
.SECONDARY:

bvp: bvp.o $(OBJECTS)
	$(FC) -o bvp bvp.o $(OBJECTS) $(LFLAGS) $(IFLAGS)

%.o: %.f90
	$(FC) -c $< $(IFLAGS)

clean:
	rm -f *.o *.mod *.png
