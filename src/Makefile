CXX=g++

FLAGS= -lz

all: main

bfgs.o : bfgs.h bfgs.cpp
	$(CXX) -c bfgs.cpp $(FLAGS)

analysisfunctions.o: analysisfunctions.h analysisfunctions.cpp
	$(CXX) -c analysisfunctions.cpp $(FLAGS)

filehandlingfunctions.o: filehandlingfunctions.h filehandlingfunctions.cpp
	$(CXX) -c filehandlingfunctions.cpp $(FLAGS)

localngsrelate-beta.o: localngsrelate-beta.cpp
	$(CXX) -c localngsrelate-beta.cpp $(FLAGS)

main: analysisfunctions.o localngsrelate-beta.o filehandlingfunctions.o bfgs.o
	$(CXX) localngsrelate.o analysisfunctions.o filehandlingfunctions.o bfgs.o  $(FLAGS) -o localngsrelate-beta

clean: 
	rm -f *.o localngsrelate-beta
