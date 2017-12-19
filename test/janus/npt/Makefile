CXX = mpicxx
CXXFLAGS = -std=c++11 -O2
CXXFLAGS += -I../../../src -I../../../../vector/src
CXXFLAGS += -DFMRES_MPI_PARALLEL
#CXXFLAGS += -DFMRES_DEBUG_PRINT

all:	main.cpp janus.hpp Makefile
	$(CXX) $(CXXFLAGS) main.cpp -o main.out
clean:
	rm main.out
