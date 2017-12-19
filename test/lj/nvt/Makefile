#CXX =g++
CXX = mpic++
CXXFLAGS = -std=c++11 -O2 -Wall
CXXFLAGS += -I../../..//src -I../../../../vector/src
CXXFLAGS += -DFMRES_MPI_PARALLEL
CXXFLAGS += -DFMRES_DEBUG_PRINT

all:	main.cpp Makefile
	$(CXX) $(CXXFLAGS) main.cpp -o main.out
