CFLAGS = -Wall -Wextra -ggdb -I/usr/include
LIBS = -L/usr/lib -lgsl -lgslcblas -lm
CC = gcc

popgen.o: popgen.c popgen.h
	$(CC) -o popgen.o -c popgen.c $(LIBS) $(CFLAGS) $(OPTS)

all: export OPTS = -O3

debug: export OPTS = -g

fastdebug: export OPTS = -g -O3

profiling: export OPTS = -g -pg

all debug fastdebug profiling: popgen.o
	$(CC) -o drifttest drifttest.c popgen.o $(LIBS) $(CFLAGS) $(OPTS)

clean:
	rm -rf *.o drifttest

