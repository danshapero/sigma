solvers ?=
linalg ?= ../
root ?= ../../../

include $(root)makefile.in

INCLUDE ?= $(root)include

solver_targets = $(solvers)cg.o $(solvers)cgs.o

$(solver_targets): $(linalg)matrix.o $(linalg)solver.o

