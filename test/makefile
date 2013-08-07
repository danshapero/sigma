test ?=
root ?= ../
src ?= ../src/

include $(root)makefile.in
include $(src)makefile

.PHONY: test

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

test_targets = $(test)type_tests $(test)util_tests $(test)graph_tests $(test)matrix_tests $(test)wrapper_tests

%_tests: libs %_tests.o
	$(FC) -o $@ $*_tests.o -lfempack -L$(LIB)
$(test)wrapper_tests: libs $(test)wrapper_tests.o
	$(CC) -o $(test)wrapper_tests $(test)wrapper_tests.o -lfempack -L$(LIB) -lgfortran
test: $(test_targets)
	$(test)./type_tests
	$(test)./util_tests
	$(test)./graph_tests
	$(test)./matrix_tests
	$(test)./wrapper_tests

