matrix ?=
graph ?= ../graph/
src ?= ../
root ?= ../../

include $(root)makefile.in
include $(graph)makefile

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

matrix_dependencies = $(src)types.o $(src)util.o $(graph_objects)
matrix_objects = $(matrix)sparse_matrices.o $(matrix)block_sparse_matrices.o $(matrix)cs_matrices.o $(matrix)coo_matrices.o $(matrix)ll_matrices.o $(matrix)ellpack_matrices.o $(matrix)bcs_matrices.o

$(matrix_objects): $(matrix_dependencies)

