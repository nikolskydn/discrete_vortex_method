CXX=g++
SRCS=$(wildcard *.cc)
HDRS=$(SRCS:%.cc=%.h) 
DLIB=libdeu2d.so

all: $(DLIB)

$(DLIB): $(SRCS) $(HDRS)
	$(CXX) -fpic -shared $(SRCS) -o $@
clean:
	rm -f *.so

