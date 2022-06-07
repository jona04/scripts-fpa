default: all

#CC = gcc -lm -W -Wall -pedantic -Wmissing-prototypes \
#	-Wredundant-decls -Wnested-externs -Wdisabled-optimization \
#	-ansi -g -fexceptions -Wno-parentheses -Wshadow -Wcast-align \
#	-Winline -Wstrict-prototypes -Wno-unknown-pragmas

# for speed:
CC = gcc -lm -O3

all: libBLG.a driver

AR = ar cr

driver: Makefile driver.c BLG.o
	$(CC) driver.c -o driver BLG.o

BLG.o: BLG.c
	$(CC) -c $*.c

libBLG.a: Makefile BLG.o
	$(AR) libBLG.a BLG.o

clean:
	rm *.o *.a driver

purge:
	rm *.o *.a driver
