# Changelog

## Version 0.3

 * New features for segmentation: penalize or ignore `None` bins
 * Segmentation accepts non-integer (i.e. normalize) count tables
 * `R/norm.R` can scale and black-list a count table
 * `R/makeNorm.R` can generate a normalization file based on other Strand-seq data (EXPERIMENTAL)
 * `R/chrom.R` highlights `None` regions 

## Version 0.2

 * Simulation accept a random generator seed for reproducible results and can output phased reads.
 * (Experimental) new HMM feature for higher ploidies
 * New plot script `R/chrom.R` to plot counts for a single chromosome, incl. SV calls and segments
 * (Experimental) new strand state classifier, because we had issues with the previously used R script.
 * Fixed bugs, e.g. in simulations (281ea3789dd2674323671bd6ab11cc118e90d7e6, abd51c38a46b3ee8c5f29849ff6daaaf245c2f80),
   in the plot function, and in some of the new functionailies

## Version 0.1

First tagged version. The tool is far from being finished, but some functionality is already available.

 * Unified interface in samtools style
   
   ```
   mosaic count
   mosaic simulate
   mosaic segment
   mosaic makebins
   ```

#### Todo: fixes

 * Don't error for wrong chromosomes during segmentation. Just skip them
 * Improve & clarify segmentation
 * Better segmentation plots (`plot_segments.R`)
 * better **NB** distribution for simulation: `boost` allows only
   integer values, `std` gives wrong numbers, but ther version in 
   `R` works corrctly! See 
   [implementation](https://github.com/wch/r-source).
 * Fine-tune HMM parameters

#### Todo: code clean up

 * Remove dependencies of HTSlib from most headers
 * More documentation
 * Fix warnings (mostly int vs unsigned comparisons)
 * Add more user input checks such as `in_range`.

#### Todo: features

 * Simulation should read + use chromosome list
 * Simulate translocations
 * Simulate LOH
 * Script to randomly generate simulation config file
 * Speed up gzip reading (maybe boost qt)
 * Count SNVs in `count` mode.
 * Add automatic tests (especially for HMM and segmentation)
