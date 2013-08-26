solver ?=
matrix ?= ../matrix/
src ?= ../
root ?= ../../

include $(root)makefile.in
include $(matrix)makefile

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

solver_dependencies = $(src)types.o $(matrix_objects)
solver_objects = $(solver)iterative_solvers.o $(solver)cg_solvers.o $(solver)bicgstab_solvers.o $(solver)jacobi_preconditioners.o

$(solver_objects): $(solver_dependencies)

