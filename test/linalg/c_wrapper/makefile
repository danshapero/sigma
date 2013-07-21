linalg_c_wrapper_tests ?=
root ?= ../../../
linalg ?= $(root)src/linalg/

include $(root)makefile.in

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

matrix_tests: $(linalg)matrix.o $(linalg_c_wrapper_tests)matrix_tests.o
	$(CC) -o $(linalg_c_wrapper_tests)matrix_tests $(linalg_c_wrapper_tests)matrix_tests.o -llinalg -lblas -llapack -lgfortran -L$(LIB)

