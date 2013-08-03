test ?=
root ?= ../
src ?= ../src/
graph_tests ?= $(test)graph/

dirs += $(graph_tests)

include $(root)makefile.in
include $(graph_tests)makefile

.PHONY: test

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

test_dependencies = $(src)types.o $(src)util.o
test_objects = $(test)type_tests.o $(test)util_tests.o
test_targets = $(test)type_tests $(test)util_tests

$(test_objects): $(test_dependencies)
$(test_targets): $(test_objects)
	$(FC) -o $(test)type_tests $(src)types.o $(test)type_tests.o -I$(INCLUDE)
	$(FC) -o $(test)util_tests $(src)util.o $(test)util_tests.o -I$(INCLUDE)
test: $(test_targets)
	$(test)./type_tests
	$(test)./util_tests
