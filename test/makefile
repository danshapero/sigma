test ?=
root ?= ../
src ?= ../src/
graph ?= $(src)graph/
matrix ?= $(src)matrix/

include $(root)makefile.in
include $(graph)makefile

.PHONY: test

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

test_dependencies = $(src)types.o $(src)util.o $(graph_objects) $(src_objects) $(matrix_objects) 
test_objects = $(test)type_tests.o $(test)util_tests.o
test_targets = $(test)type_tests $(test)util_tests $(test)graph_tests $(test)matrix_tests

$(test_objects): $(test_dependencies)
$(test)type_tests: $(src)types.o $(test)type_tests.o
	$(FC) -o $(test)type_tests $(test)type_tests.o $(src)types.o
$(test)util_tests: $(src)util.o $(test)util_tests.o
	$(FC) -o $(test)util_tests $(test)util_tests.o $(src)util.o
$(test)graph_tests: $(test)graph_tests.o $(test_dependencies)
	$(FC) -o $(test)graph_tests $(test)graph_tests.o $(test_dependencies)
$(test)matrix_tests: $(test)matrix_tests.o $(test_dependencies)
	$(FC) -o $(test)matrix_tests $(test)matrix_tests.o $(test_dependencies)
test: $(graph_tests_dependencies) $(test_objects) $(test_targets)
	$(test)./type_tests
	$(test)./util_tests
	$(test)./graph_tests
	$(test)./matrix_tests

