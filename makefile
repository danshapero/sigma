include makefile.in

.SECONDARY:
.PHONY: clean

INCLUDE = include/
LIB = lib/

root =
src = src/
test = test/
dirs = $(src) $(test)

include $(src)makefile
include $(test)makefile


clean:
	@for dir in $(dirs); \
	do \
		rm -f $$dir$ *.o ; \
	done; \
	rm -f $(INCLUDE)*.mod
