solvers ?=
include $(solvers)../../../makefile.in

F90INCLUDE ?= $(solvers)../../../include

solver_targets = $(solvers)cg.o $(solvers)cgs.o

$(solver_targets): $(solvers)../matrix.o $(solvers)../solver.o

