STATIC ?= 0
DEBUG ?= 0


# Submodules
PWD = $(shell pwd)
HTSLIB_ROOT ?= ${PWD}/src/htslib/


# Flags
CXX = g++
CXXFLAGS += -isystem ${HTSLIB_ROOT} -std=c++11 -pedantic -W -Wall -Wno-unknown-pragmas -D__STDC_LIMIT_MACROS
LDFLAGS += -L${HTSLIB_ROOT}


# Additional flags for release/debug
ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz
else
	LDFLAGS += -lhts -lz -Wl,-rpath,${HTSLIB_ROOT}
endif
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	# DNDEBUG removes the macro "assert" completely; not used in my code
	CXXFLAGS += -O3 -DNDEBUG 
endif


# Sources
MAINSOURCES = src/main.cpp $(wildcard src/*.h)
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)


# Targets
TARGETS = .htslib src/main

all: $(TARGETS)

src/main: $(MAINSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

.htslib: $(HTSLIBSOURCES)
	cd src/htslib && make && make lib-static && cd ../../ && touch .htslib

clean:
	cd src/htslib && make clean
	rm -f $(TARGETS) $(TARGETS:=.o)