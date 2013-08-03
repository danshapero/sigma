test ?=
root ?= ../
src ?= ../src/
graph_tests ?= $(test)graph/

dirs += $(graph_tests)

include $(root)makefile.in
#include $(graph_tests)makefile

.PHONY: test

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

test_dependencies = $(src)types.o $(src)util.o
graph_tests_dependencies = $(test_dependencies) $(graph_objects)
test_objects = $(test)type_tests.o $(test)util_tests.o
graph_tests_objects = $(test)graph_tests.o
test_targets = $(test)type_tests $(test)util_tests $(test)graph_tests

$(test_objects): $(test_dependencies)
$(graph_tests_objects): $(graph_tests_dependencies)
$(test_targets): $(test_objects) $(graph_tests_objects)
	$(FC) -o $(test)type_tests $(src)types.o $(test)type_tests.o -I$(INCLUDE)
	$(FC) -o $(test)util_tests $(src)util.o $(test)util_tests.o -I$(INCLUDE)
	$(FC) -o $(test)graph_tests $(test)graph_tests.o $(graph_tests_dependencies) -I$(INCLUDE)
test: $(test_targets)
	$(test)./type_tests
	$(test)./util_tests
	$(test)./graph_tests

