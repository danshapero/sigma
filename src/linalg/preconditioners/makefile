preconditioners ?=
linalg ?= ../
matrix_formats ?= $(linalg)matrix_formats/
root ?= ../../../

include $(root)makefile.in

INCLUDE ?= $(root)include

preconditioner_targets = $(preconditioners)jacobi.o $(preconditioners)bjacobi.o $(preconditioners)ilu.o

$(preconditioner_targets): $(linalg)matrix.o $(linalg)solver.o

ilu.o: $(linalg)matrix.o $(linalg)solver.o $(matrix_formats)csr.o