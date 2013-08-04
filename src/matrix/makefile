matrix ?=
graph ?= ../graph/
src ?= ../
root ?= ../../

include $(root)makefile.in

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

matrix_dependencies = $(src)util.o $(graph_objects)
matrix_objects = $(matrix)sparse_matrices.o $(matrix)csr_matrices.o

$(matrix_objects): $(matrix_dependencies)

