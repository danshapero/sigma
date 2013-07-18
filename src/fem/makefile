fem ?=
mesh ?= ../mesh/
linalg ?= ../linalg/

include $(fem)../../makefile.in

INCLUDE ?= $(fem)../../include/
LIB ?= $(fem)../../lib/

libfem.a: $(mesh)meshes.o $(linalg)matrix.o $(fem)fem.o
	$(AR) rc $(LIB)libfem.a $(fem)fem.o
	$(RANLIB) $(LIB)libfem.a

