graph_tests ?=
root ?= ../../
src ?= $(root)src/
graph ?= $(src)graph/

.PHONY: graph_tests

include $(root)makefile.in
include $(graph)makefile

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

graph_tests_dependencies = $(src)types.o $(src)util.o $(graph_objects)
graph_tests_objects = $(graph_tests)coo_graph_tests.o $(graph_tests)ll_graph_tests.o
graph_tests_targets = $(graph_tests)coo_graph_tests $(graph_test)ll_graph_tests

$(graph_tests_objects): $(graph_tests_dependencies)
$(graph_tests_targets): $(graph_tests_objects)
	$(FC) -o $(graph_tests)coo_graph_tests $(graph_tests)coo_graph_tests.o $(graph)coo_graphs.o $(graph)abstract_graphs.o -I$(INCLUDE)
	$(FC) -o $(graph_tests)ll_graph_tests $(graph_tests)ll_graph_tests.o $(graph)ll_graphs.o $(graph)abstract_graphs.o $(src)types.o -I$(INCLUDE)
graph_tests: $(graph_tests_targets)
	$(graph_tests)./coo_graph_tests
	$(graph_tests)./ll_graph_tests

