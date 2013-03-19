OBJECTS = abstract_matrix_mod.o csr_matrix_mod.o linalg_mod.o
DFLAGS = -g -fbounds-check

FC = gfortran

.PHONY: clean, all
.SECONDARY:

all: $(OBJECTS)

%.o: %.f90
	$(FC) -c $<

clean:
	rm -f *.o *.mod
