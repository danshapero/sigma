OBJECTS = mesh_mod.o
DFLAGS = -g -fbounds-check -Wall

FC = gfortran

.PHONY:
.SECONDARY:

lib: $(OBJECTS)
	ar cru libmesh.a $(OBJECTS)
	ranlib libmesh.a
	mv libmesh.a ../../lib
	mv *.mod ../../include

%.o: %.f90
	$(FC) -c -fPIC -Wall $<

clean:
	rm -f *.o *.mod
