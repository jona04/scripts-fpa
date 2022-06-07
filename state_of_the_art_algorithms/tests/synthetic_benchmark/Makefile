# Simple makefile to compile the synthetic (random) tests

# Compilation flags
CC=gcc
FC=gfortran
OPTFLAGS=-march=native -O3 -ffast-math
CFLAGS=-ansi -std=c99 -Wall -pedantic $(OPTFLAGS)
ifeq ($(CFG),debug)
	CFLAGS += -DDEBUG
endif
FFLAGS=$(OPTFLAGS)

# Standard libraries
LIBS=-lm

########## If you do not have BLAS available comment below

# System BLAS linking information
CFLAGS += -DLAPACK
LIBS += -lblas

########## end BLAS section

# Targets
BINARY=synthetic
OBJS=fixed_point.o synthetic_benchmark.o cont_quad_knapsack.o third_party_methods.o dsel05.o dSFMT.o
all: $(BINARY)

# Dependencies and build rules
synthetic_benchmark.o: synthetic_benchmark.c ../../lib/cont_quad_knapsack.h ../../lib/third_party_methods.h ../../lib/fixed_point.h 
	$(CC) $(CFLAGS) -DDSFMT_MEXP=19937 -c synthetic_benchmark.c

cont_quad_knapsack.o: ../../lib/cont_quad_knapsack.h ../../lib/cont_quad_knapsack.c
	$(CC) $(CFLAGS) -c ../../lib/cont_quad_knapsack.c

third_party_methods.o: ../../lib/third_party_methods.h ../../lib/third_party_methods.c ../../lib/cont_quad_knapsack.h
	$(CC) $(CFLAGS) -c ../../lib/third_party_methods.c

fixed_point.o: ../../lib/fixed_point.h ../../lib/fixed_point.c ../../lib/cont_quad_knapsack.h
	$(CC) $(CFLAGS) -c ../../lib/fixed_point.c

dsel05.o: ../../../third_party/dsel05.f
	$(FC) $(FFLAGS) -c ../../../third_party/dsel05.f

dSFMT.o: ../../../third_party/dSFMT-src-2.2.2/dSFMT.c
	$(CC) -O3 -fno-strict-aliasing -DDSFMT_MEXP=19937 -c ../../../third_party/dSFMT-src-2.2.2/dSFMT.c

synthetic: $(OBJS)
	$(CC) $(CFLAGS) -o synthetic $(OBJS) $(LIBS)

clean:
	rm $(OBJS) $(BINARY) *~ times.dat
