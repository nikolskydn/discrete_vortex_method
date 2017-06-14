CXX=g++
CXXFLAGS=-std=c++11 -fopenmp
SRCS=$(wildcard make*.cc)
TRGS=$(SRCS:%.cc=%)
LDIR=../lib
FIGSRC=figure2d.cc
all: $(TRGS)


$(TRGS): %: %.cc 
	make -C $(LDIR)
	$(CXX) $(CXXFLAGS) -o $@ $(FIGSRC) $< -L$(LDIR) -ldeu2d

clean:
	rm -f $(TRGS)
