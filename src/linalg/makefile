linalg ?=
root ?= ../../

include $(root)makefile.in

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

matrix_formats = $(linalg)matrix_formats/
include $(matrix_formats)makefile

solvers = $(linalg)solvers/
include $(solvers)makefile

preconditioners = $(linalg)preconditioners/
include $(preconditioners)makefile

c_wrapper = $(linalg)c_wrapper/
include $(c_wrapper)makefile

dirs += $(matrix_formats) $(solvers) $(preconditioners) $(c_wrapper)

linalg_targets = $(linalg)matrix.o $(linalg)solver.o $(linalg)permutation.o $(matrix_format_targets) $(solver_targets) $(preconditioner_targets) $(c_wrapper_targets)

linalg.o: $(linalg_targets)

liblinalg.a: $(linalg_targets) $(linalg)linalg.o
	$(AR) rc $(LIB)liblinalg.a $(linalg)linalg.o $(linalg_targets)
	$(RANLIB) $(LIB)liblinalg.a

