#
# Copyright (C) 2017 Sascha Meiers
# Distributed under the MIT software license, see the accompanying
# file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.

library(data.table)
library(assertthat)
library(ggplot2)
library(scales)
library(cowplot)

args <- commandArgs(trailingOnly = T)

if (length(args) !=4 | !grepl('\\.pdf$', args[4])) {
    warning("Usage: Rscript R/plot_segments.R count-file segment-file chrom output-pdf")
    quit(status = 1)
}
f_in <- args[1]
f_seg <- args[2]
chrom_ <- args[3]
pdf_out <- args[4]

format_Mb <- function(x) {
    paste(comma(x/1e6), "Mb")
}

# if gzip
if (substr(f_in,nchar(f_in)-2,nchar(f_in)) == ".gz")
    f_in = paste("zcat",f_in)

# strip "chr" from chromosomes
chrom_ = sub('^chr', '', chrom_)

# Read counts & filter chromosomes (this is human-specific)
d = fread(f_in)

# Check that correct files are given:
assert_that("chrom" %in% colnames(d))
assert_that("start" %in% colnames(d) && is.integer(d$start))
assert_that("end" %in% colnames(d) && is.integer(d$end))
assert_that("sample" %in% colnames(d))
assert_that("cell" %in% colnames(d))
assert_that("w" %in% colnames(d) && is.integer(d$w))
assert_that("c" %in% colnames(d) && is.integer(d$c))

# Re-name and -order chromosomes
d = d[, chrom := sub('^chr','',chrom)][]


# subset to chrom_
assert_that(chrom_ %in% d$chrom)
d = d[chrom == chrom_,]


# Read segmentation
seg = fread(f_seg)
seg = seg[, chrom := sub('^chr','',chrom)][]
seg = seg[chrom == chrom_,]
# distribute colors in seg
seg <- seg[, color := as.factor( (0:(k-1) %% 3)), by = k][]


# Subset to fewer cells
max_cells = 12
N_cells = nrow(unique(d[, .(sample,cell)]))
if (N_cells > max_cells) {
    message("Subsetting to first ", max_cells, " cells")
    d <- d[paste(sample,cell) %in% unique(d[, paste(sample,cell)])[1:max_cells], ]
    N_cells = max_cells
}


cairo_pdf(pdf_out, width=10, height=2+N_cells/2, onefile = T)
for(k_ in unique(seg$k))
{
    message("Plotting ", N_cells, " cells with ", k_, " breakpoints")
    plt <- ggplot(d) + 
        geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax = -w), fill='sandybrown') +
        geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax =  c), fill='paleturquoise4') +
        facet_grid(sample+cell~., scales = "free_y") +
        ylab("Watson | Crick") + xlab(NULL) +
        scale_x_continuous(breaks = pretty_breaks(12), labels = format_Mb) +
        scale_y_continuous(breaks = pretty_breaks(3)) + 
        theme(panel.spacing = unit(0.2, "lines"),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.background = element_rect(fill = NA, colour=NA)) + 
        guides(fill = FALSE) +
        geom_rect(data = seg[k == k_,], aes(xmin = start, xmax=end, ymin=-Inf, ymax = Inf, fill = color), alpha = 0.33) +
        scale_fill_brewer(type='seq') +
        ggtitle(paste("Chromosome", chrom_, "with", k_, "breakpoints"))

    print(plt)
}
dev.off()
