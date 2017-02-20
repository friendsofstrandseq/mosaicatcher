STATIC ?= 0
DEBUG ?= 0


# Submodules
PWD = $(shell pwd)
HTSLIB_ROOT ?= ${PWD}/src/htslib/
BOOST_ROOT ?= ${PWD}/src/boost/
STOCHHMM_ROOT ?= ${PWD}/src/StochHMM/


# Flags
CXX = g++
CXXFLAGS += -fdiagnostics-color -isystem ${HTSLIB_ROOT} -isystem ${STOCHHMM_ROOT} -isystem ${BOOST_ROOT} -std=c++11 -pedantic -W -Wall -Wno-unknown-pragmas -D__STDC_LIMIT_MACROS
LDFLAGS += -L${HTSLIB_ROOT} -L${BOOST_ROOT}/stage/lib -L${STOCHHMM_ROOT}/src -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time


# Additional flags for release/debug
ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz
else
	LDFLAGS += -lhts -lz -Wl,-rpath,${HTSLIB_ROOT},-rpath,${BOOST_ROOT}/stage/lib
endif
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	# DNDEBUG removes the macro "assert" completely; not used in my code
	# Using O2 instead of O3 because of a segfault. Maybe relted to https://github.com/samtools/htslib/issues/400
	CXXFLAGS += -O2 -DNDEBUG -fno-tree-vectorize
endif


# Sources
MAINSOURCES = src/main.cpp $(wildcard src/*.h)
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
BOOSTSOURCES = $(wildcard src/boost/libs/iostreams/include/boost/iostreams/*.hpp)
STOCHHMMSOURCES = $(wildcard src/StochHMM/)


# Targets
TARGETS = .htslib .boost .stochHMM src/main src/test_hmm

all: $(TARGETS)

src/main: $(MAINSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/test_hmm: .stochHMM src/test_hmm.cpp
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS) ${STOCHHMM_ROOT}/src/libstochhmm.a

.htslib: $(HTSLIBSOURCES)
	cd src/htslib && make && make lib-static && cd ../../ && touch .htslib

.boost: $(BOOSTSOURCES)
	cd src/boost && ./bootstrap.sh --prefix=${PWD}/src/boost --without-icu --with-libraries=iostreams,filesystem,system,program_options,date_time && ./b2 && ./b2 headers && cd ../../ && touch .boost

.stochHMM: $(STOCHHMMSOURCES)
	cd src/StochHMM && ./configure && make && cd ../../ && touch .stochHMM

clean:
	cd src/htslib && make clean
	cd src/boost && ./b2 --clean-all
	cd src/stochHMM && make clean
	rm -f $(TARGETS) $(TARGETS:=.o)

