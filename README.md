# Processing Strand-seq data

## Installation (beta)

Mosaicatcher can be built using [Cmake](https://cmake.org/) (v3.0) on **Linux** and **MacOS**. 

It relies on two external dependecies

 * [boost libraries](http://www.boost.org/) >= 1.50. This needs to be installed on your system
 * [HTSlib](https://github.com/samtools/htslib) >= 1.3.1. Cmake should be able to install this for you

```
git clone https://github.com/friendsofstrandseq/mosaicatcher.git --branch develop
cd mosaicatcher
mkdir build
cd build
cmake ../src
make
./mosaic --version
```

## Strand-seq read counting and generation of QC plots

Mosaicatcher counts Strand-seq reads and classifies strand states of each chromosome in each cell
using a hidden Markov model.

Choose between bins of fixed width (`-w`) or predefined bins (`-b`).

```
./src/main -o counts.txt.gz -i counts.info -x data/exclude/GRCh38_full_analysis_set_plus_decoy_hla.exclude -w 200000 cell1.bam cell2.bam [...]
```

To generate QC plots from these tables run

```
Rscript R/qc.R counts.txt.gz counts.info counts.pdf
```

## Data input

Sequencing reads from a single cell should be grouped into a single BAM file and all reads within that BAM file will be treated as one single cell.

The BAM file must contain a single read group. Cells are grouped to samples based on this `SM` tag.


## Simulation

Simulate strand-seq data and SVs on the level of binned counts.

```
./src/simul -o data/simulation/out.txt.gz data/simulation/example.txt
Rscript R/qc.R data/simulation/out.txt.gz data/simulation/out.txt.pdf
```


