# Simple makefile to compile the synthetic (random) tests

# Compilation flags
CC=gcc-12
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
OBJS=fixed_point.o fixed_point_new.o fixed_point_a.o fixed_point_ai.o synthetic_benchmark.o cont_quad_knapsack.o third_party_methods.o dsel05.o dSFMT.o
all: $(BINARY)

# Dependencies and build rules
synthetic_benchmark.o: synthetic_benchmark.c ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h ../../../state_of_the_art_algorithms/lib/third_party_methods.h ../../../state_of_the_art_algorithms/lib/fixed_point.h ../../../state_of_the_art_algorithms/lib/fixed_point_new.h ../../lib/fixed_point_a.h ../../lib/fixed_point_ai.h 
	$(CC) $(CFLAGS) -DDSFMT_MEXP=19937 -c synthetic_benchmark.c

cont_quad_knapsack.o: ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.c
	$(CC) $(CFLAGS) -c ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.c

third_party_methods.o: ../../../state_of_the_art_algorithms/lib/third_party_methods.h ../../../state_of_the_art_algorithms/lib/third_party_methods.c ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h
	$(CC) $(CFLAGS) -c ../../../state_of_the_art_algorithms/lib/third_party_methods.c

fixed_point.o: ../../../state_of_the_art_algorithms/lib/fixed_point.h ../../../state_of_the_art_algorithms/lib/fixed_point.c ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h
	$(CC) $(CFLAGS) -c ../../../state_of_the_art_algorithms/lib/fixed_point.c

fixed_point_new.o: ../../../state_of_the_art_algorithms/lib/fixed_point_new.h ../../../state_of_the_art_algorithms/lib/fixed_point_new.c ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h
	$(CC) $(CFLAGS) -c ../../../state_of_the_art_algorithms/lib/fixed_point_new.c

fixed_point_a.o: ../../lib/fixed_point_a.h ../../lib/fixed_point_a.c ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h
	$(CC) $(CFLAGS) -c ../../lib/fixed_point_a.c

fixed_point_ai.o: ../../lib/fixed_point_ai.h ../../lib/fixed_point_ai.c ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h
	$(CC) $(CFLAGS) -c ../../lib/fixed_point_ai.c

dsel05.o: ../../../third_party/dsel05.f
	$(FC) $(FFLAGS) -c ../../../third_party/dsel05.f

dSFMT.o: ../../../third_party/dSFMT-src-2.2.2/dSFMT.c
	$(CC) -O3 -fno-strict-aliasing -DDSFMT_MEXP=19937 -c ../../../third_party/dSFMT-src-2.2.2/dSFMT.c

synthetic: $(OBJS)
	$(CC) $(CFLAGS) -o synthetic $(OBJS) $(LIBS)

clean:
	rm $(OBJS) $(BINARY) *~ times.dat
