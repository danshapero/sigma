preconditioners ?=
include $(preconditioners)../../../makefile.in

F90INCLUDE ?= $(preconditioners)../../../include

preconditioner_targets = $(preconditioners)jacobi.o $(preconditioners)bjacobi.o $(preconditioners)ilu.o

$(preconditioner_targets): $(preconditioners)../matrix.o $(preconditioners)../solver.o

