CC = gcc
LD = gcc
CFLAGS = -O3 -fopenmp   
LDFLAGS = -fopenmp -llapacke -lm 

OBJFILES = time.o
MAIN = time

all: $(MAIN)

$(MAIN): $(OBJFILE)
	$(LD) -o $(OBJFILE) $(LDFLAGS)

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) 

clean:
	rm -f $(OBJFILES) $(SHARED)
