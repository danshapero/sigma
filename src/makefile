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
src_objects = $(src)fempack.o

$(src_objects): $(src_dependencies)

