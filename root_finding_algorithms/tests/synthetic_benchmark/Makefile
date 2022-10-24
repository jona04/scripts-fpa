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
OBJS=fixed_point_new.o synthetic_benchmark.o bissection_method.o regula_falsi.o secant_method.o dsel05.o dSFMT.o
all: $(BINARY)

# Dependencies and build rules
synthetic_benchmark.o: synthetic_benchmark.c ../../lib/bissection_method.h ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h 
	$(CC) $(CFLAGS) -DDSFMT_MEXP=19937 -c synthetic_benchmark.c

bissection_method.o: ../../lib/bissection_method.h ../../lib/bissection_method.c ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h 
	$(CC) $(CFLAGS) -c ../../lib/bissection_method.c

regula_falsi.o: ../../lib/bissection_method.h ../../lib/regula_falsi.c ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h 
	$(CC) $(CFLAGS) -c ../../lib/regula_falsi.c

fixed_point_new.o: ../../../state_of_the_art_algorithms/lib/fixed_point_new.h ../../../state_of_the_art_algorithms/lib/fixed_point_new.c ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h
	$(CC) $(CFLAGS) -c ../../../state_of_the_art_algorithms/lib/fixed_point_new.c

secant_method.o: ../../lib/bissection_method.h ../../lib/secant_method.c ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h
	$(CC) $(CFLAGS) -c ../../lib/secant_method.c

dsel05.o: ../../../third_party/dsel05.f
	$(FC) $(FFLAGS) -c ../../../third_party/dsel05.f

dSFMT.o: ../../../third_party/dSFMT-src-2.2.2/dSFMT.c
	$(CC) -O3 -fno-strict-aliasing -DDSFMT_MEXP=19937 -c ../../../third_party/dSFMT-src-2.2.2/dSFMT.c

synthetic: $(OBJS)
	$(CC) $(CFLAGS) -o synthetic $(OBJS) $(LIBS)

clean:
	rm $(OBJS) $(BINARY) *~ times.dat

# # Simple makefile to compile the synthetic (random) tests

# # Compilation flags
# CC=gcc-12
# FC=gfortran
# OPTFLAGS=-march=native -O3 -ffast-math
# CFLAGS=-ansi -std=c99 -Wall -pedantic $(OPTFLAGS)
# ifeq ($(CFG),debug)
# 	CFLAGS += -DDEBUG
# endif
# FFLAGS=$(OPTFLAGS)

# # Standard libraries
# LIBS=-lm

# ########## If you do not have BLAS available comment below

# # System BLAS linking information
# CFLAGS += -DLAPACK
# LIBS += -lblas

# ########## end BLAS section

# # Targets
# BINARY=synthetic
# OBJS=fixed_point_new.o synthetic_benchmark.o bissection_method.o regula_falsi.o secant_method.o dsel05.o dSFMT.o
# all: $(BINARY)

# # Dependencies and build rules
# synthetic_benchmark.o: synthetic_benchmark.c ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h ../../../state_of_the_art_algorithms/lib/third_party_methods.h ../../lib/bissection_method.h ../../lib/regula_falsi.h ../../lib/secant_method.h
# 	$(CC) $(CFLAGS) -DDSFMT_MEXP=19937 -c synthetic_benchmark.c

# bissection_method.o: ../../lib/bissection_method.h ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h ../../lib/bissection_method.c
# 	$(CC) $(CFLAGS) -c ../../lib/bissection_method.c

# regula_falsi.o: ../../lib/bissection_method.h ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h ../../lib/regula_falsi.c
# 	$(CC) $(CFLAGS) -c ../../lib/regula_falsi.c

# fixed_point_new.o: ../../../state_of_the_art_algorithms/lib/fixed_point_new.h ../../../state_of_the_art_algorithms/lib/fixed_point_new.c ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h
# 	$(CC) $(CFLAGS) -c ../../../state_of_the_art_algorithms/lib/fixed_point_new.c

# secant_method.o: ../../lib/bissection_method.h ../../../state_of_the_art_algorithms/lib/cont_quad_knapsack.h ../../lib/secant_method.c
# 	$(CC) $(CFLAGS) -c ../../lib/secant_method.c

# dsel05.o: ../../../third_party/dsel05.f
# 	$(FC) $(FFLAGS) -c ../../../third_party/dsel05.f

# dSFMT.o: ../../../third_party/dSFMT-src-2.2.2/dSFMT.c
# 	$(CC) -O3 -fno-strict-aliasing -DDSFMT_MEXP=19937 -c ../../../third_party/dSFMT-src-2.2.2/dSFMT.c

# synthetic: $(OBJS)
# 	$(CC) $(CFLAGS) -o synthetic $(OBJS) $(LIBS)

# clean:
# 	rm $(OBJS) $(BINARY) *~ times.dat
