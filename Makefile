CXX = g++
CXXFLAGS += -std=c++11 -pedantic -W -Wall -Wno-unknown-pragmas -D__STDC_LIMIT_MACROS 
LDFLAGS += 


# Sources
MAINSOURCES = src/main.cpp $(wildcard src/*.h)

# Targets
TARGETS = src/main

all: $(TARGETS)

src/main: $(MAINSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGETS) $(TARGETS:=.o)