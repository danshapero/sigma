OBJECTS = fem_mod.o
IFLAGS = -I../../include
DFLAGS = -g -fbounds-check -Wall

FC = gfortran

.PHONY: clean
.SECONDARY:

lib: $(OBJECTS)
	ar cru libfem.a $(OBJECTS)
	ranlib libfem.a
	mv libfem.a ../../lib
	mv *.mod ../../include

%.o: %.f90
	$(FC) -c -fPIC -Wall $< $(IFLAGS)

clean:
	rm -f *.o *.mod
