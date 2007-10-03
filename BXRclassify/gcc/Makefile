CFLAG=-DUSE_GCC -g -I/home/lshenzhi/BXRclassify/src/tclap-1.0.5/include 

all: BXRclassify
BXRclassify: MultiModel.o BXRclassify.o parseModelFiles.o readLabelFile.o readTestVectorFile.o genResultFile.o
	g++ $(CFLAG) -static -o BXRclassify MultiModel.o BXRclassify.o parseModelFiles.o  readLabelFile.o readTestVectorFile.o genResultFile.o
	cp BXRclassify ../src/

BXRclassify.o: ../src/BXRclassify.cpp
	g++ $(CFLAG) -c $<

MultiModel.o: ../src/MultiModel.cpp ../src/MultiModel.h 
	g++ $(CFLAG) -c ../src/MultiModel.cpp -o MultiModel.o 

parseModelFiles.o: ../src/parseModelFiles.cpp ../src/MultiModel.h
	g++ $(CFLAG) -c ../src/parseModelFiles.cpp -o parseModelFiles.o

readLabelFile.o: ../src/readLabelFile.cpp ../src/MultiModel.h
	g++ $(CFLAG) -c ../src/readLabelFile.cpp -o readLabelFile.o

readTestVectorFile.o: ../src/readTestVectorFile.cpp ../src/MultiModel.h
	g++ $(CFLAG) -c ../src/readTestVectorFile.cpp -o readTestVectorFile.o

genResultFile.o: ../src/genResultFile.cpp ../src/MultiModel.h
	g++ $(CFLAG) -c ../src/genResultFile.cpp -o genResultFile.o


clean:
	rm *.o BXRclassify ../src/BXRclassify
