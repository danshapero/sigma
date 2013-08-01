matrix ?=
graph ?= ../graph/
src ?= ../
root ?= ../../

include $(root)makefile.in

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

matrix_dependencies = $(src)util.o $(graph_objects)
matrix_objects = $(matrix)abstract_sparse_matrices.o

$(matrix_objects): $(matrix_dependencies)

