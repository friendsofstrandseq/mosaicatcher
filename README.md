# Processing Strand-seq data

Structural Variant calling from single-cell Strand-seq* data.

This software is in development.

**Falconer E et al., 2012 ([PMID 23042453](https://www.ncbi.nlm.nih.gov/pubmed/23042453))*


## Installation

Mosaicatcher can be built using [Cmake](https://cmake.org/) (v3.0) on **Linux** and **MacOS**. 

It relies on two external dependecies

 * [boost libraries](http://www.boost.org/) >= 1.50. This needs to be installed on your system
 * [HTSlib](https://github.com/samtools/htslib) >= 1.3.1. Cmake should be able to install this for you

```
git clone https://github.com/friendsofstrandseq/mosaicatcher.git
cd mosaicatcher
mkdir build
cd build
cmake ../src
make
./mosaic --version
```

## Strand-seq read counting and plotting

Mosaicatcher counts Strand-seq reads and classifies strand states of each chromosome in each cell
using a Hidden Markov Model.

Choose between bins of fixed width (`-w`) or predefined bins (`-b`). 
Here is an example for bins with a fixed width of 200kb:

```
./build/mosaic count \
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

Sequencing reads should be supplied in exactly one BAM file per single cell.
Each BAM file must contain a single read group. Cells are grouped to 
*samples* based on their `SM` tag.


## Strand-seq simulations

Simulate strand-seq data and SVs on the level of binned counts. You are asked to specify an *SV config* file such as in the example `data/simulation/example.txt`.

Then run

```
./build/mosaic simulate \
    -o counts.txt.gz \
	 svconfig.txt
Rscript R/qc.R counts.txt.gz counts.pdf
```


