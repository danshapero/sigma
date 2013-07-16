include makefile.in

.SECONDARY:
.PHONY: clean veryclean

INCLUDE = include/
LIB = lib/

dirs = src/linalg/ src/mesh/ src/fem/

fem = src/fem/
include $(fem)makefile

linalg = src/linalg/
include $(linalg)makefile

mesh = src/mesh/
include $(mesh)makefile

examples = examples/
include $(examples)makefile

libs: liblinalg.a libmesh.a libfem.a

clean:
	@for dir in $(dirs); \
	do \
		rm -f $$dir$ *.o ; \
	done

veryclean:
	@for dir in $(dirs); \
	do \
		rm -f $$dir$ *.o ; \
	done; \
	rm -f $(LIB)*.a $(INCLUDE)*.mod

