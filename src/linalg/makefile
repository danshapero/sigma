linalg ?=
include $(linalg)../../makefile.in

INCLUDE ?= $(linalg)../../include/
LIB ?= $(linalg)../../lib/

dirs += $(linalg)matrix_formats/ $(linalg)solvers/ $(linalg)preconditioners/

matrix_formats = $(linalg)matrix_formats/
include $(matrix_formats)makefile

solvers = $(linalg)solvers/
include $(solvers)makefile

preconditioners = $(linalg)preconditioners/
include $(preconditioners)makefile

linalg_targets = $(linalg)matrix.o $(linalg)solver.o $(linalg)permutation.o $(matrix_format_targets) $(solver_targets) $(preconditioner_targets)

linalg.o: $(linalg_targets)

liblinalg.a: $(linalg_targets) $(linalg)linalg.o
	$(AR) rc $(LIB)liblinalg.a $(linalg)linalg.o $(linalg_targets)
	$(RANLIB) $(LIB)liblinalg.a

