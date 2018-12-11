#===General Variables===#
CC=gcc
CFLAGS=-Wall -Wextra -g3

all: makeAll

makeAll: makeHelper makeRing makeStats makeRandom makeRandom makeGenetic makeMain
	$(CC) $(CFLAGS) helper.o ring.o stats.o random.o genetic.o main.o -o genetic -ldl -lm -lblas -llapack -lpthread 

makeMain: main.c
	$(CC) $(CFLAGS) -c main.c -o main.o 

makeGenetic: genetic.c genetic.h
	$(CC) $(CFLAGS) -c genetic.c -o genetic.o 

makeRandom: random_numbers.c random_numbers.h
	$(CC) $(CFLAGS) -c random_numbers.c -o random.o 

makeStats: stats.c stats.h
	$(CC) $(CFLAGS) -c stats.c -o stats.o 

makeRing: ring.c ring.h
	$(CC) $(CFLAGS) -c ring.c -o ring.o

makeHelper: helper.c helper.h macros.h
	$(CC) $(CFLAGS) -c helper.c -o helper.o


.PHONY: clean

clean:
	rm -f *~ $(ODIR)/*.o $(IDIR)/*~
