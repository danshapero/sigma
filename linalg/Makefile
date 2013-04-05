OBJECTS = sparse_matrix_mod.o csr_matrix_mod.o iterative_solver_mod.o cg_solver_mod.o jacobi_mod.o
DFLAGS = -g -fbounds-check -Wall

FC = gfortran

.PHONY: clean, all
.SECONDARY:

lib: $(OBJECTS)
	ar cru liblinalg.a $(OBJECTS)
	ranlib liblinalg.a

%.o: %.f90
	$(FC) -c -fPIC -Wall $<

clean:
	rm -f *.o *.mod
