# Changelog


## Version 0.2

Small features:
 * A seed to the simulation random generator for reproducible results
 * (Experimental) new HMM feature for higher ploidies

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

 * `git push --delete segmentation`
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
