##############################################################################
### Compiler

CC = gcc

### Compiler flags

CFLAGS  = -Ofast -march=native -fopenmp
#CFLAGS  = -g -Wall -fopenmp

### Linked libraries

LDLIBS = -lmitmojco -fopenmp -lm

### make all

.PHONY: all
all: example-1 example-2 example-3 example-4 example-5

DEPENDS = Makefile

### make example-1

example-1: example-1.o
	$(CC) $@.o -o $@ $(LDLIBS)

example-1.o: $(DEPENDS)

### make example-2

example-2: example-2.o
	$(CC) $@.o -o $@ $(LDLIBS)

example-2.o: $(DEPENDS)

### make example-3

example-3: example-3.o
	$(CC) $@.o -o $@ $(LDLIBS)

example-3.o: $(DEPENDS)

### make example-4

example-4: example-4.o
	$(CC) $@.o -o $@ $(LDLIBS)

example-4.o: $(DEPENDS)

### make example-5

example-5: example-5.o
	$(CC) $@.o -o $@ $(LDLIBS)

example-5.o: $(DEPENDS)

### make clean

.PHONY: clean
clean:
	rm -f example-1 example-1.o
	rm -f example-2 example-2.o
	rm -f example-3 example-3.o
	rm -f example-4 example-4.o
	rm -f example-5 example-5.o

