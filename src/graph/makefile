graph ?=
src ?= ../
root ?= ../../

include $(root)makefile.in

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

graph_targets = $(graph)abstract_graphs.o $(graph)coo_graphs.o $(graph)cs_graphs.o

graph_targets: $(src)utility_algorithms.o $(graph_targets)

