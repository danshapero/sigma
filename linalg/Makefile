OBJECTS = sparse_matrix_mod.o csr_matrix_mod.o iterative_solver_mod.o cg_solver_mod.o
DFLAGS = -g -fbounds-check

FC = gfortran

.PHONY: clean, all
.SECONDARY:

all: $(OBJECTS)

%.o: %.f90
	$(FC) -c $<

clean:
	rm -f *.o *.mod
