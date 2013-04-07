OBJECTS = sparse_matrix_mod.o csr_matrix_mod.o iterative_solver_mod.o cg_solver_mod.o jacobi_mod.o nullpc_mod.o
DFLAGS = -g -fbounds-check -Wall

FC = gfortran

.PHONY: clean
.SECONDARY:

lib: $(OBJECTS)
	ar cru liblinalg.a $(OBJECTS)
	ranlib liblinalg.a
	mv liblinalg.a ../../lib
	mv *.mod ../../include
	rm *.o

%.o: %.f90
	$(FC) -c -fPIC $<

clean:
	rm -f *.o *.mod
