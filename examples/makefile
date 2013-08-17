root ?= ../
examples ?=
src ?= ../src/

include $(root)makefile.in
include $(src)makefile

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

poisson: libs fem.o poisson.o
	$(FC) -o poisson poisson.o fem.o -lfempack -llapack -L$(LIB)

