# Makefile for Assignment 2, DTU 02635, Fall 2017
CC=gcc
CFLAGS=-Wall -std=c99 -g
EXECUTABLE=lssolve

# Detect OS and set linker flags accordingly
ifdef SystemRoot
 	# Windows
	LDFLAGS=-L.
	LDLIBS=-lopenblas
else
 	ifeq ($(shell uname), Darwin)
		# macOS
		LDLIBS=-llapack -lblas
	else
		# Gbar
		LDFLAGS=-L/usr/lib64/atlas
		LDLIBS=-llapack -lf77blas -latlas
	endif
endif

.PHONY: all clean run

all: $(EXECUTABLE)

$(EXECUTABLE): matrix_io.o lssolve.o

lssolve.o: lssolve.h
matrix_io.o: matrix_io.h

clean:
	-$(RM) *.o $(EXECUTABLE)

run: $(EXECUTABLE)
	./$< $(ARGS)
