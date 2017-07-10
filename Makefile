##############################################################################
### Compiler

CC = gcc

### Compiler flags

CFLAGS  = -Ofast -march=native -fopenmp
#CFLAGS  = -g -Wall -fopenmp

### Linked libraries

LDLIBS = -lm -fopenmp

### make all

.PHONY: all
all: example-1 example-2 example-3 example-4 example-5

mitmojco.o: mitmojco.h Makefile
opt_filter.o: opt_filter.h Makefile
DEPENDS = mitmojco.h mitmojco.c opt_filter.h opt_filter.c Makefile

### make example-1

example-1: example-1.o mitmojco.o opt_filter.o
	$(CC) $@.o mitmojco.o opt_filter.o -o $@ $(LDLIBS)

example-1.o: $(DEPENDS)

### make example-2

example-2: example-2.o mitmojco.o opt_filter.o
	$(CC) $@.o mitmojco.o opt_filter.o -o $@ $(LDLIBS)

example-2.o: $(DEPENDS)

### make example-3

example-3: example-3.o mitmojco.o opt_filter.o
	$(CC) $@.o mitmojco.o opt_filter.o -o $@ $(LDLIBS)

example-3.o: $(DEPENDS)

### make example-4

example-4: example-4.o mitmojco.o opt_filter.o
	$(CC) $@.o mitmojco.o opt_filter.o -o $@ $(LDLIBS)

example-4.o: $(DEPENDS)

### make example-5

example-5: example-5.o mitmojco.o opt_filter.o
	$(CC) $@.o mitmojco.o opt_filter.o -o $@ $(LDLIBS)

example-5.o: $(DEPENDS)

### make clean

.PHONY: clean
clean:
	rm -f mitmojco.o opt_filter.o
	rm -f example-1 example-1.o
	rm -f example-2 example-2.o
	rm -f example-3 example-3.o
	rm -f example-4 example-4.o
	rm -f example-5 example-5.o

