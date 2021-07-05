CC = gcc
CC2 = g++
LD = gcc
LD2 = g++
CFLAGS = -O2 -fopenmp -fPIC 
LDFLAGS = -fopenmp -llapacke -llapack -lblas -lm -lquadmath
all = time
all: $(all)

time: time.o
	$(LD) -o time time.o $(LDFLAGS) -pg

time.o: time.c
	$(CC) -c -o time.o time.c $(CFLAGS) -pg

clean:
	rm -f *.o time
