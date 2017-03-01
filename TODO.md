## How to compile on EMBL server:

```
ssh pc-korbel17

```

# To do

## Build/make:

  * put `#define BOOST_DISABLE_ASSERTS` before boost inclusions
  * give header files a namespace and global macro against double-inclusion
  * remove `-fdiagnostics-color` flag from CXXFLAGS
  * DNDEBUG removes the macro "assert" completely; not used in my code
  * Enable `-O3` once it runs without segfault. Might be related to https://github.com/samtools/htslib/issues/400

## Optimisations:

  * unused chromosomes and blacklisted bins could be skipped
  
  
