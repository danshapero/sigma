graph ?=
src ?= ../
root ?= ../../

include $(root)makefile.in

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

graph_dependencies = $(src)types.o $(src)util.o
graph_objects = $(graph)abstract_graphs.o $(graph)coo_graphs.o $(graph)cs_graphs.o

$(graph_objects): $(graph_dependencies)

