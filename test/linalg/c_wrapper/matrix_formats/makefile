c_matrix_format_tests ?= 
root ?= ../../../../
linalg ?= $(root)src/linalg/
matrix_formats ?= $(linalg)matrix_formats/
c_matrix_formats ?= $(linalg)c_wrapper/matrix_formats/

include $(root)makefile.in

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

csr_tests: $(linalg)matrix.o $(matrix_formats)csr.o $(c_matrix_formats)csr_c.o $(c_matrix_format_tests)csr_tests.o
	$(CC) -o $(c_matrix_format_tests)csr_tests $(c_matrix_formats)csr_c.o $(c_matrix_format_tests)csr_tests.o -lgfortran -llinalg -L$(LIB)

