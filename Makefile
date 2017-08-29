STATIC ?= 0
DEBUG ?= 0


# Submodules
PWD = $(shell pwd)
HTSLIB_ROOT ?= ${PWD}/src/htslib/
BOOST_ROOT ?= ${PWD}/src/boost/


# Flags
CXX ?= g++
CXXFLAGS += -isystem ${HTSLIB_ROOT} -isystem ${BOOST_ROOT} -std=c++11 -pedantic -W -Wall -Wno-unknown-pragmas -D__STDC_LIMIT_MACROS
LDFLAGS += -L${HTSLIB_ROOT} -L${BOOST_ROOT}/stage/lib -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time


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
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
BOOSTSOURCES = $(wildcard src/boost/libs/iostreams/include/boost/iostreams/*.hpp)


# Targets
TARGETS = .htslib .boost src/main src/calc_bins src/segmentation

all: $(TARGETS)

src/main: .boost .htslib src/main.cpp $(wildcard src/*.hpp)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/calc_bins: .boost .htslib src/calc_bins.cpp
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/segmentation: .boost .htslib src/segmentation.cpp
        $(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

.htslib: $(HTSLIBSOURCES)
	cd src/htslib && make && make lib-static && cd ../../ && touch .htslib

.boost: $(BOOSTSOURCES)
	cd src/boost && ./bootstrap.sh --prefix=${PWD}/src/boost --without-icu --with-libraries=iostreams,filesystem,system,program_options,date_time && ./b2 && ./b2 headers && cd ../../ && touch .boost

doc/html/index.html: doc/Doxyfile
	cd src && doxygen ../doc/Doxyfile

clean:
	cd src/htslib && make clean
	cd src/boost && ./b2 --clean-all
	rm -f $(TARGETS) $(TARGETS:=.o)

