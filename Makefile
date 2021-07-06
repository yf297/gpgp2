CC = gcc
LD = gcc
CFLAGS = -O2 -fopenmp -fPIC 
LDFLAGS = -fopenmp -llapacke -llapack -lblas -lm -lquadmath
all = time
all: $(all)

time: time.o
	$(LD) -o time time.o $(LDFLAGS) 

time.o: time.c
	$(CC) -c -o time.o time.c $(CFLAGS)

clean:
	rm -f *.o time
