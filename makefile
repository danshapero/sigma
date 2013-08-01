include makefile.in

.SECONDARY:
.PHONY: clean

INCLUDE = include/
LIB = lib/

root =
src = src/
include $(src)makefile

test = test/
include $(test)makefile

dirs = $(src) $(test)


clean:
	@for dir in $(dirs); \
	do \
		rm -f $$dir$ *.o ; \
	done; \
	rm -f $(INCLUDE)*.mod
