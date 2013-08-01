matrix ?=
graph ?= ../graph/
src ?= ../
root ?= ../../

include $(root)makefile.in

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

matrix_targets = $(matrix)abstract_sparse_matrices.o

matrix_targets: $(src)utility_algorithms.o $(graph_targets) $(matrix_targets)

