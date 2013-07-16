mesh ?=
include $(mesh)../../makefile.in

INCLUDE ?= $(mesh)../../include/
LIB ?= $(mesh)../../lib/

libmesh.a: $(mesh)meshes.o
	$(AR) rc $(LIB)libmesh.a $(mesh)meshes.o
	$(RANLIB) $(LIB)libmesh.a

