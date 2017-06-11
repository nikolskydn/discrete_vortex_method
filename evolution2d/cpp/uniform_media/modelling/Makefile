CXX=g++
CXXFLAGS=-std=c++11 -fopenmp
SRCS=$(wildcard *.cc)
TRGS=$(SRCS:%.cc=%)
LDIR=../lib
all: $(TRGS)

$(TRGS): %: %.cc
	make -C $(LDIR)
	$(CXX) $(CXXFLAGS) -o $@ $< -L$(LDIR) -ldeu2d

clean:
	rm -f $(TRGS)

