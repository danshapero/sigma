OBJECTS = mesh_mod.o linalg/abstract_matrix_mod.o linalg/csr_matrix_mod.o fem_mod.o linalg/linalg_mod.o #permutation_mod.o
LFLAGS = -L/home/daniel/programs/netcdf/lib -lnetcdf -lnetcdff
IFLAGS = -I/home/daniel/programs/netcdf/include -I./linalg
DFLAGS = -g -fbounds-check

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

%.exe: $(OBJECTS) %.o
#ifeq ($(WITH_DEBUG),0)
	$(FC) $*.o -o $@ $(OBJECTS) $(LFLAGS) $(IFLAGS)
#else
#	$(FC) $*.o -o $@ $(OBJECTS) $(LFLAGS) $(IFLAGS) $(DFLAGS)
#endif

linalg/%.o: linalg/%.f90
	cd linalg && make all

%.o: %.f90
#ifeq ($(WITH_DEBUG),0)
	$(FC) -c $< $(IFLAGS)
#else
#	$(FC) -c $< $(IFLAGS) $(DFLAGS)
#endif

clean:
	rm -f *.o *.mod *.exe *.png
