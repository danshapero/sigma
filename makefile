include makefile.in

.SECONDARY:
.PHONY: clean

INCLUDE = include/
LIB = lib/

dirs = $(src)

root =
src = src/
include $(src)makefile


clean:
	@for dir in $(dirs); \
	do \
		rm -f $$dir$ *.o ; \
	done; \
	rm -f $(INCLUDE)*.mod
