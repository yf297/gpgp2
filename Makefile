CC = clang
CFLAGS = -fopenmp -fpic -O3 -Wall
LDFLAGS = -shared -fopenmp 
LIBS = -lm
LIBBLAS = `pkg-config --libs openblas`
LIBLAPACK = `pkg-config --libs lapack`
LIBLAPACKE = `pkg-config --libs lapacke`
INCBLAS = `pkg-config --cflags openblas`

OBJFILES = vecchia.o
SHARED = vecchia

all: $(SHARED)

$(SHARED): $(OBJFILES)
	$(CC) $(LDFLAGS) $(LIBS) $(LIBLAPACKE) $(LIBLAPACK) $(LIBBLAS) -o $(SHARED) $(OBJFILES) 

%.o: %.c
	$(CC) $(CFLAGS) $(INCBLAS) -c -o $@ $<

clean:
	rm -f $(OBJFILES) $(SHARED)
