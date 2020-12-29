CC = gcc
LD = gcc
CFLAGS = -fopenmp -O3
LDFLAGS = -fopenmp -shared  -llapacke 

OBJFILES = vecchia.o
SHARED = vecchia

all: $(SHARED)

$(SHARED): $(OBJFILES)
	$(LD) -o $(SHARED) $(OBJFILES) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJFILES) $(SHARED)
