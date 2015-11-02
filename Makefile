CC=g++
LIBDIR=/home/sgeraedt/Loops/
CFLAGS=-c -O3 -fopenmp -I$(LIBDIR) 
LDFLAGS=-fopenmp -I$(LIBDIR) 
SOURCES=lattice.cpp lattice.h ising.cpp ising.h mainLoopy.cpp loopy.cpp loopy.h VTable.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=a.out

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm *.o
