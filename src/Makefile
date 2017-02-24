CFLAGS = -Wall -Wextra -ggdb -I/usr/include
LIBS = -L/usr/lib -lgsl -lgslcblas -lm
CC = gcc

popgen.o: popgen.c popgen.h
	$(CC) -o popgen.o -c popgen.c $(LIBS) $(CFLAGS) $(OPTS)

readdata.o: readdata.c readdata.h
	$(CC) -o readdata.o -c readdata.c $(LIBS) $(CFLAGS) $(OPTS)

test.o: test.c test.h 
	$(CC) -o test.o -c test.c  $(LIBS) $(CFLAGS) $(OPTS)

all: export OPTS = -O3

debug: export OPTS = -g

fastdebug: export OPTS = -g -O3

profiling: export OPTS = -g -pg

all debug fastdebug profiling: popgen.o readdata.o test.o
	$(CC) -o ../bin/drifttest main.c popgen.o readdata.o test.o $(LIBS) $(CFLAGS) $(OPTS)

clean:
	rm -rf *.o ../bin/drifttest

