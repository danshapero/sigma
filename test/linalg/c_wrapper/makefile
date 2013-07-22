c_wrapper_tests ?=
root ?= ../../../
linalg ?= $(root)src/linalg/

include $(root)makefile.in

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

matrix_tests: $(linalg)matrix.o $(c_wrapper_tests)matrix_tests.o
	$(CC) -o $(c_wrapper_tests)matrix_tests $(c_wrapper_tests)matrix_tests.o -llinalg -lblas -llapack -lgfortran -L$(LIB)

solver_tests: $(linalg)solver.o $(linalg)c_wrapper/solver_c.o $(c_wrapper_tests)solver_tests.o
	$(CC) -o $(c_wrapper_tests)solver_tests $(c_wrapper_tests)solver_tests.o -llinalg -lblas -llapack -lgfortran -L$(LIB)

