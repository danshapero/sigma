src ?=
root ?= ../

include $(root)makefile.in

INCLUDE ?= $(root)include/
LIB ?= $(root)lib/

graph ?= $(src)graph/
matrix ?= $(src)matrix/

dirs += $(graph) $(matrix)

include $(graph)makefile
include $(matrix)/makefile

src_dependencies = $(graph_objects) $(matrix_objects)
src_objects = $(src)wrapper.o $(src)fempack.o

libs: $(src)util.o $(src)types.o $(src_dependencies) $(src_objects)
	$(AR) rc $(LIB)libfempack.a $(src)util.o $(src)types.o $(graph_objects) $(matrix_objects) $(src)wrapper.o $(src)fempack.o
	$(RANLIB) $(LIB)libfempack.a
