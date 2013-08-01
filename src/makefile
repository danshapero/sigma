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

src_targets = $(src)graphs.o

src_targets: $(graph_targets) $(src)graphs.o $(matrix_targets)

