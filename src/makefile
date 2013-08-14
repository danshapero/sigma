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
src_objects = $(src)conversions.o $(src)wrapper.o $(src)meshes.o $(src)fempack.o

libs: $(src)types.o $(src)util.o $(src_dependencies) $(src_objects)
	$(AR) rc $(LIB)libfempack.a $(src)util.o $(src)types.o $(graph_objects) $(matrix_objects) $(src_objects)
	$(RANLIB) $(LIB)libfempack.a
