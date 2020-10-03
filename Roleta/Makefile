# Name of the project
PROJ_NAME=roulette_kcluster

# Flags CPLEX
# CPLEXFLAGS=-m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD -I/opt/ibm/ILOG/CPLEX_Studio129/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio129/concert/include  -L/opt/ibm/ILOG/CPLEX_Studio129/cplex/lib/x86-64_linux/static_pic -lilocplex -lcplex -L/opt/ibm/ILOG/CPLEX_Studio129/concert/lib/x86-64_linux/static_pic -lconcert -lm -pthread -std=c++0x -ldl

# Flags CPP
CFLAGS=-std=c++11 -static-libstdc++ -static-libgcc -Wall -O3

# Source files
DECODERS_SOURCE=$(wildcard ../Decoders/*.cpp)
TIMER_SOURCE=$(wildcard ../Timer/*.cpp)
IN_SOURCE=$(wildcard *.cpp)

ALL_SOURCE=$(DECODERS_SOURCE) $(TIMER_SOURCE) $(IN_SOURCE) 

# Object files
OBJ=$(ALL_SOURCE:.cpp=.o)

all:$(OBJ)
	g++  $(OBJ) -o $(PROJ_NAME) $(CPLEXFLAGS) $(CFLAGS)

%.o: %.cpp %.hpp
	g++ -c $< -o $@ $(CFLAGS)

main.o: main.cpp
	g++ -c -o main.o main.cpp $(CPLEXFLAGS) $(CFLAGS)

clean:
	rm -rf ../Decoders/*.o
	rm -rf ../Timer/*.o
	rm -rf *.o $(PROJ_NAME) *~
