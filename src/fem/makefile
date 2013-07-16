fem ?= 
include $(fem)../../makefile.in

INCLUDE ?= $(fem)../../include/
LIB ?= $(fem)../../lib/

mesh = $(fem)../mesh/
#include $(mesh)makefile

linalg = $(fem)../linalg/
#include $(linalg)makefile

libfem.a: $(mesh)meshes.o $(linalg)matrix.o $(fem)fem.o
	$(AR) rc $(LIB)libfem.a $(fem)fem.o
	$(RANLIB) $(LIB)libfem.a

