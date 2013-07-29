matrix_format_tests ?=
root ?= ../../../
linalg ?= $(root)src/linalg/
matrix_formats ?= $(linalg)matrix_formats/

include $(root)makefile.in

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

csr_tests: $(linalg)matrix.o $(matrix_format_tests)csr_tests.o
	$(FC) -o $(matrix_format_tests)csr_tests $(matrix_format_tests)csr_tests.o -llinalg -lblas -llapack -L$(LIB)

bsr_tests: $(linalg)matrix.o $(matrix_format_tests)bsr_tests.o
	$(FC) -o $(matrix_format_tests)bsr_tests $(matrix_format_tests)bsr_tests.o -llinalg -lblas -llapack -L$(LIB)

ellpack_tests: $(linalg)matrix.o $(matrix_format_tests)ellpack_tests.o
	$(FC) -o $(matrix_format_tests)ellpack_tests $(matrix_format_tests)ellpack_tests.o -llinalg -lblas -llapack -L$(LIB)

