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
OBJS=fixed_point.o synthetic_benchmark.o bissection_method.o regula_falsi.o secant_method.o dsel05.o dSFMT.o
all: $(BINARY)

# Dependencies and build rules
synthetic_benchmark.o: synthetic_benchmark.c ../../lib/bissection_method.h ../../lib/fixed_point.h 
	$(CC) $(CFLAGS) -DDSFMT_MEXP=19937 -c synthetic_benchmark.c

bissection_method.o: ../../lib/bissection_method.h ../../lib/fixed_point.h ../../lib/bissection_method.c
	$(CC) $(CFLAGS) -c ../../lib/bissection_method.c

regula_falsi.o: ../../lib/bissection_method.h ../../lib/fixed_point.h ../../lib/regula_falsi.c
	$(CC) $(CFLAGS) -c ../../lib/regula_falsi.c

fixed_point.o: ../../lib/fixed_point.h ../../lib/fixed_point.c
	$(CC) $(CFLAGS) -c ../../lib/fixed_point.c

secant_method.o: ../../lib/bissection_method.h ../../lib/fixed_point.h ../../lib/secant_method.c
	$(CC) $(CFLAGS) -c ../../lib/secant_method.c

dsel05.o: ../../../third_party/dsel05.f
	$(FC) $(FFLAGS) -c ../../../third_party/dsel05.f

dSFMT.o: ../../../third_party/dSFMT-src-2.2.2/dSFMT.c
	$(CC) -O3 -fno-strict-aliasing -DDSFMT_MEXP=19937 -c ../../../third_party/dSFMT-src-2.2.2/dSFMT.c

synthetic: $(OBJS)
	$(CC) $(CFLAGS) -o synthetic $(OBJS) $(LIBS)

clean:
	rm $(OBJS) $(BINARY) *~ times.dat
