OBJECTS = mesh_mod.o fem_mod.o #linalg/sparse_matrix_mod.o linalg/csr_matrix_mod.o fem_mod.o linalg/iterative_solver_mod.o linalg/cg_solver_mod.o linalg/jacobi_mod.o #permutation_mod.o
LFLAGS = -lnetcdf -lnetcdff -L./linalg -llinalg
IFLAGS = -I/usr/local/include -I./linalg
DFLAGS = -g -fbounds-check -Wall

MESH = circle.1
RHS = none
BND = none
SOL = u
PIC = $(SOL)
MODE = robin
WITH_PLOTS = 1
#WITH_DEBUG = 0

FC = gfortran

.PHONY: clean
.SECONDARY:

poisson: poisson.exe
	./poisson.exe meshes/$(MESH) $(SOL) $(RHS) $(BND) $(MODE)
ifeq ($(WITH_PLOTS),1)
	python plotting.py $(MESH) $(SOL).nc $(PIC)
endif

%.exe: linalg/liblinalg.a $(OBJECTS) %.o
ifeq ($(WITH_DEBUG),0)
	$(FC) $*.o -o $@ $(OBJECTS) $(LFLAGS) $(IFLAGS)
else
	$(FC) $*.o -o $@ $(OBJECTS) $(LFLAGS) $(IFLAGS) $(DFLAGS)
endif

linalg/liblinalg.a: linalg/sparse_matrix_mod.f90
	cd linalg && make lib

%.o: %.f90
ifeq ($(WITH_DEBUG),0)
	$(FC) -c $< $(IFLAGS)
else
	$(FC) -c $< $(IFLAGS) $(DFLAGS)
endif

clean:
	rm -f *.o *.mod *.exe *.png
