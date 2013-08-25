graph ?=
src ?= ../
root ?= ../../

include $(root)makefile.in

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

graph_dependencies = $(src)types.o $(src)util.o
graph_objects = $(graph)graphs.o $(graph)coo_graphs.o $(graph)cs_graphs.o $(graph)ellpack_graphs.o $(graph)ll_graphs.o

$(graph_objects): $(graph_dependencies)

