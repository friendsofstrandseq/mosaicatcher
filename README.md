[![Build Status](https://travis-ci.org/friendsofstrandseq/mosaicatcher.svg?branch=develop)](https://travis-ci.org/friendsofstrandseq/mosaicatcher)

# Processing Strand-seq data

This software is part of a larger [pipeline](https://github.com/friendsofstrandseq/pipeline) to call structural variants in single-cell Strand-seq data.

> For optimal integration with [pipeline](https://github.com/friendsofstrandseq/pipeline) version 1.0, please use version `0.3.1-dev`.


## Installation

Mosaicatcher can be built using [Cmake](https://cmake.org/) (v3.0) on **Linux** and **MacOS**. 

It relies on two external dependecies, which should both be installed on your system:

 * [boost libraries](http://www.boost.org/) >= 1.50.
 * [HTSlib](https://github.com/samtools/htslib) >= 1.3.1.

```
git clone https://github.com/friendsofstrandseq/mosaicatcher.git
cd mosaicatcher
mkdir build
cd build
cmake ../src
make
./mosaicatcher --version
```

## Strand-seq read counting and plotting

Mosaicatcher counts Strand-seq reads and classifies strand states of each chromosome in each cell
using a Hidden Markov Model.

Choose between bins of fixed width (`-w`) or predefined bins (`-b`). 
Here is an example for bins with a fixed width of 200kb:

```
./build/mosaicatcher count \
    -o counts.txt.gz \
    -i counts.info \
    -x data/exclude/GRCh38_full_analysis_set_plus_decoy_hla.exclude \
    -w 200000 \
    cell1.bam cell2.bam [...]
```

To generate QC plots from these tables run

```
Rscript R/qc.R \
    counts.txt.gz \
    counts.info \
    counts.pdf
```

### Data input

* Sequencing reads should be supplied in exactly one BAM file per single cell
* Each BAM file must contain a single read group (`@RG`). Cells are grouped into *samples* by using the same `SM` tag.
* BAM files must be sorted and indexed.


## Strand-seq simulations

Simulate strand-seq data and SVs on the level of binned counts. You are asked to specify an *SV config* file such as in the example `data/simulation/example.txt`.

Then run

```
./build/mosaicatcher simulate \
    -o counts.txt.gz \
	 svconfig.txt
Rscript R/qc.R counts.txt.gz counts.pdf
```


## References

For information on Strand-seq see

> Falconer E *et al.*, 2012 (doi: [10.1038/nmeth.2206](https://doi.org/10.1038/nmeth.2206))
