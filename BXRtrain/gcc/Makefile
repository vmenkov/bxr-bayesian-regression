# C++ objects

CC=g++
CFLAGS=-O0 -c -g -Wno-deprecated -I/home/lshenzhi/BXR/BXRtrain/src/tclap-1.0.5/include -I../src -DUSE_GCC -DSYMBOLIC
LDFLAGS=-G
SOURCES=DataFactory.cpp Matrix.cpp StrataSplit.cpp Milliseconds.cpp OutputFormat.cpp Likelihood.cpp PolyZO.cpp QNOpt.cpp CommandLine.cpp
CSOURCES=blmvm-1.1ez.c
OBJECTS=blmvm-1.1ez.o DataFactory.o Matrix.o StrataSplit.o Milliseconds.o OutputFormat.o Likelihood.o PolyZO.o QNOpt.o CommandLine.o

EXECUTABLE=BXRtrain

all: $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) -static $(OBJECTS) -o $@

%.o: ../src/%.cpp
	$(CC) $(CFLAGS) $< 


%.o: ../src/%.c
	$(CC) -x c $(CFLAGS) $< 

clean:
	rm *.o
	rm BXRtrain

tar:
	tar cvf BXRtrain.tar *.cpp *.h *.c Makefile
	gzip -9 BXRtrain.tar

